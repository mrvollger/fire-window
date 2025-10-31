#!/usr/bin/env Rscript

# load libraries
suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(ggpubr)
    library(glue)
    library(cowplot)
})

# command line args
parser <- ArgumentParser(description = "Generate actuation plots from pileup files")
parser$add_argument("BED",
    type = "character",
    help = "Path to manifest BED file.
#chrom: chromosome/contig name, start: region start position (0-based)\n
end: region end position
center: centered position for plotting along that region
pileup: path to fibertools pileup BED file (must be tabix-indexed)
"
)

parser$add_argument("-o", "--output",
    type = "character", required = TRUE,
    help = "Output PDF file path"
)
parser$add_argument("-t", "--height-per-track",
    type = "double", default = 1,
    help = "Height per track in inches [default: %(default)s]"
)
parser$add_argument("-w", "--width",
    type = "double", default = 8,
    help = "Figure width in inches [default: %(default)s]"
)
parser$add_argument("-c", "--color",
    type = "character", default = "darkred",
    help = "Color for the ribbon plot [default: %(default)s]"
)
parser$add_argument("-b", "--highlight-bed",
    type = "character", default = NULL,
    help = "Optional BED file with regions to highlight (chrom, start, end)"
)
parser$add_argument("--highlight-color",
    type = "character", default = "darkorange",
    help = "Color for highlighted regions [default: %(default)s]"
)

args <- parser$parse_args()


Height_per_track <- args$height_per_track
Width <- args$width
Manifest <- args$BED
Output <- args$output
Color <- args$color
Highlight_bed <- args$highlight_bed
Highlight_color <- args$highlight_color

manifest_df <- fread(Manifest) %>%
    mutate(
        chrom = `#chrom`,
        region = paste0(chrom, ":", start + 1, "-", end),
        region_length = end - start,
        region_start = start,
        region_end = end,
        region_center = center
    ) %>%
    select(
        -start, -end, -center
    )
# if strand is not in df set it to +
if (!"strand" %in% colnames(manifest_df)) {
    manifest_df$strand <- "+"
}

read_pileup_region <- function(pileup_file, region) {
    cmd <- glue("tabix -h {pileup_file} {region}")
    fread(cmd = cmd, sep = "\t", header = TRUE) %>%
        mutate(
            chrom = `#chrom`,
            region = region
        )
}

# for every pileup in manifest_df read in the pileup file for the region
all_pileups <- bind_rows(lapply(1:nrow(manifest_df), function(i) {
    row <- manifest_df[i, ]
    df <- read_pileup_region(row$pileup, row$region) %>%
        mutate(
            sample_id = row$sample_id,
        )
    df
})) %>%
    merge(
        manifest_df %>%
            select(-chrom, -`#chrom`),
        by = "region",
        all.x = TRUE
    ) %>%
    mutate(
        tstart = start,
        tend = end,
        start = ifelse(strand == "+", tstart, region_length - tend),
        end = ifelse(strand == "+", tend, region_length - tstart),
        # Strand correction for region boundaries and center
        temp_region_start = region_start,
        temp_region_end = region_end,
        region_start = ifelse(strand == "+", temp_region_start, region_length - temp_region_end),
        region_end = ifelse(strand == "+", temp_region_end, region_length - temp_region_start),
        region_center = ifelse(strand == "+", region_center, region_length - region_center)
    ) %>%
    select(-tstart, -tend, -temp_region_start, -temp_region_end) %>%
    mutate(
        region = factor(region, levels = rev(manifest_df$region)),
        frac = fire_coverage / coverage,
    )

# center the pileups around the center
center_pileups <- all_pileups %>%
    mutate(
        start = start - region_center,
        end = end - region_center,
        region_start = region_start - region_center,
        region_end = region_end - region_center,
    )

# read and process highlight regions if provided
highlight_regions <- NULL
if (!is.null(Highlight_bed)) {
    highlight_regions <- fread(Highlight_bed, select = 1:3) %>%
        setNames(c("chrom", "start", "end")) %>%
        merge(
            manifest_df %>%
                select(chrom, region, region_length, strand,
                       region_start, region_end, region_center),
            by = "chrom",
            allow.cartesian = TRUE
        ) %>%
        # filter to only highlights that overlap the region
        filter(
            start < region_end & end > region_start
        ) %>%
        mutate(
            # clip to region boundaries
            start = pmax(start, region_start),
            end = pmin(end, region_end),
            # apply strand correction
            temp_start = start,
            temp_end = end,
            start = ifelse(strand == "+", temp_start,
                          region_length - temp_end),
            end = ifelse(strand == "+", temp_end,
                        region_length - temp_start),
            # center around region_center
            start = start - region_center,
            end = end - region_center
        ) %>%
        select(region, start, end) %>%
        mutate(
            region = factor(region, levels = levels(center_pileups$region))
        )
}

# make a plot of the actuation faceting on the region
center_pileups %>%
    # to plot this correctly I need to copy each row twice
    # once for start and once for end
    bind_rows(
        center_pileups %>%
            mutate(
                end = start
            )
    ) %>%
    ggplot() +
    geom_ribbon(
        aes(
            x = start,
            ymin = 0,
            ymax = frac,
        ),
        # alpha=0.2
        fill = Color,
    ) +
    # add highlighted regions if provided
    {
        if (!is.null(highlight_regions)) {
            geom_rect(
                data = highlight_regions,
                aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
                fill = Highlight_color,
                alpha = 0.3,
                inherit.aes = FALSE
            )
        }
    } +
    facet_wrap(region ~ ., ncol = 1) +
    # draw a vertical line at the centering point
    geom_vline(
        data = . %>%
            select(region) %>%
            distinct(),
        aes(xintercept = 0),
        color = "black",
        linetype = "dashed",
    ) +
    # draw a horizonal line for the region length
    geom_segment(
        data = . %>%
            select(region, region_start, region_end) %>%
            distinct(),
        aes(
            x = region_start,
            xend = region_end,
            y = 0,
            yend = 0
        ),
        color = "darkblue",
        linewidth = 0.2
    ) +
    scale_y_continuous(
        "% actuation",
        limits = c(0, 1),
        labels = scales::percent,
    ) +
    scale_x_continuous(
        "Position (bp)",
        label = scales::comma,
        # expand = c(0, 0)
    ) +
    theme_cowplot(
        font_size = 8
    )

ggsave(
    Output,
    width = Width,
    height = Height_per_track * nlevels(center_pileups$region)
)
