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
parser = ArgumentParser(description="Generate actuation plots from pileup files")
parser$add_argument("BED", type="character",
                    help="Path to manifest BED file.
#chrom: chromosome/contig name, start: region start position (0-based)\n
end: region end position
center: centered position for plotting along that region
pileup: path to fibertools pileup BED file (must be tabix-indexed)
")

parser$add_argument("-o", "--output", type="character", required=TRUE,
                    help="Output PDF file path")
parser$add_argument("-t", "--height-per-track", type="double", default=1,
                    help="Height per track in inches [default: %(default)s]")
parser$add_argument("-w", "--width", type="double", default=8,
                    help="Figure width in inches [default: %(default)s]")
parser$add_argument("-c", "--color", type="character", default="darkred",
                    help="Color for the ribbon plot [default: %(default)s]")

args = parser$parse_args()




Height_per_track = args$height_per_track
Width = args$width
Manifest = args$BED
Output = args$output
Color = args$color

manifest_df = fread(Manifest) %>%
    mutate(
        chrom=`#chrom`,
        region = paste0(chrom, ":", start+1, "-", end),
        region_length = end - start,
        region_start = start,
        region_end = end,
        region_center = center
    ) %>%
    select(
        -start, -end, -center
    )
# if strand is not in df set it to +
if(!"strand" %in% colnames(manifest_df)) {
    manifest_df$strand = "+"
}

read_pileup_region = function(pileup_file, region) {
    cmd = glue("tabix -h {pileup_file} {region}")
    fread(cmd=cmd, sep="\t", header=TRUE) %>%
        mutate(
            chrom=`#chrom`,
            region = region
        )
}

# for every pileup in manifest_df read in the pileup file for the region
all_pileups = bind_rows(lapply(1:nrow(manifest_df), function(i) {
    row = manifest_df[i,]
    df = read_pileup_region(row$pileup, row$region) %>%
        mutate(
            sample_id = row$sample_id,
        )
    df
})) %>%
    merge(
        manifest_df %>%
            select(-chrom, -`#chrom`),
        by="region",
        all.x=TRUE
    ) %>%
    mutate(
        tstart = start,
        tend = end,
        start = ifelse(strand == "+", tstart, region_length - tend),
        end = ifelse(strand == "+", tend, region_length - tstart)
    ) %>%
    select(-tstart, -tend) %>%
    mutate(
        region = factor(region, levels=rev(manifest_df$region)),
        frac = fire_coverage / coverage,
    )

# center the pileups around the center
center_pileups = all_pileups %>%
    mutate(
        start = start - region_center,
        end = end - region_center,
        region_start = region_start - region_center,
        region_end = region_end - region_center,
    )

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
            x=start,
            ymin=0,
            ymax=frac,
        ),
        #alpha=0.2
        fill=Color,
    ) +
    facet_wrap(region ~ ., ncol=1) +
    # draw a vertical line at the centering point
    geom_vline(
        data= . %>%
            select(region) %>%
            distinct(),
        aes(xintercept = 0),
        color="black",
        linetype="dashed",
    ) +
    # draw a horizonal line for the region length
    geom_segment(
        data= . %>%
            select(region, region_start, region_end) %>%
            distinct(),
        aes(
            x=region_start,
            xend=region_end,
            y=0,
            yend=0
        ),
        color="darkblue",
        linewidth=0.2
    ) +
    scale_y_continuous(
        "% actuation",
        limits=c(0,1),
        labels = scales::percent,
    ) +
    scale_x_continuous(
        "Position (bp)",
        label = scales::comma,
        #expand = c(0, 0)
    ) +
    theme_cowplot(
        font_size = 8
    )

ggsave(
    Output,
    width=Width,
    height=Height_per_track * nlevels(center_pileups$region)
)
