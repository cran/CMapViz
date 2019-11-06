#' @export
#'
#' @import ggplot2
#' @import readxl
#' @import dplyr
#' @import reshape2
#' @import scales
#' @import stringr
#' @title Bubble plot of CMap output table
#' @description
#' This function allows the user to represent the Connectivity Map (CMap) result table (broadinstitute)
#' under the form of a bubble plot representing statistics and cell lines:
#' - each drug is represented along the y axis according to its enrichment value
#' - each drug is represented along the x axis according to the cell line tested and
#' within the cell line according to batch specificity (0-50% <- vertical line -> 50-100%)
#' @name bubble_plot
#' @rdname bubble_plot
#' @aliases bubble_plot
#' @usage
#' bubble_plot(path, plot, enrichment, abs.enrich.cutoff=NULL, n.rep.cutoff=NULL ,
#'         jittering=FALSE, return.gg.table= FALSE, output_path = NULL)
#' @param path path of the excel file
#' @param plot what data to plot: molecules only (plot="molecules") or molecules by cell lines batch (plot="cell.lines")
#' @param enrichment whether to plot positive or negative enrichment
#' @param abs.enrich.cutoff minimum value of enrichment to include a batch
#' @param n.rep.cutoff  minimum number of replicates to include a batch (default=NULL)
#' @param jittering whether apply jittering to the values to avoid points overlap
#' @param return.gg.table table prepared for ggplot, allows the user to customize the graphical representation
#' @param output_path path for the experiment output folder, returns data table and figure (default=NULL)
#'
#' @examples
#' file.path <- system.file("extdata", "example.xls", package = "CMapViz")
#' #display results by cell lines, with negative enrichment (absolute cutoff: 0.5), and at least n=5.
#' #molecule position with respect of dotted line is the specificity of the molecule itself:
#' #left side of dotted line if specificity < 50 or right side of dotted line if specificity > 50 )
#' bubble_plot(file.path,
#'     plot = "cell.lines", enrichment = "negative", abs.enrich.cutoff = 0.5,
#'     n.rep.cutoff = 5, output_path = NULL
#' )
#' @return ggplot object - bubble plot



bubble_plot <- function(path,
                        plot,
                        enrichment,
                        abs.enrich.cutoff = NULL,
                        n.rep.cutoff = NULL,
                        jittering = FALSE,
                        return.gg.table = FALSE,
                        output_path = NULL) {

    # import the specified data
    if (plot == "cell.lines") {
        expected=c("rank","cmap.name.and.cell.line","mean","n","enrichment",
                "p","specificity","percent non-null")
        sheet <- 2         #"by cmap name and cell line"
        table <- read_excel(path, sheet = sheet)

        #error handling sheet column names
        if (!identical(colnames(table),expected))
        {
                stop('The input file is not properly formated:
                -please use the original output excel table from the Connectivity Map (CMap)
                -second-sheet column names shoud be: \n\t\t',paste(expected,collapse=', '))
        }
        names <- table$cmap.name.and.cell.line
        names_and_cl <- colsplit(names, "-(?!.*-)", c("drug", "cell.line"))
        rank <- rep(1, dim(table)[1])
        df <- cbind(names_and_cl, table[3:dim(table)[2]], rank)
    } else if (plot == "molecules") {
        expected=c("rank","cmap name","mean","n","enrichment",
                "p","specificity","percent non-null")
        # sheet <- "by cmap name"
        sheet = 1
        table <- read_excel(path, sheet = sheet)

        #error handling sheet column names
        if (!identical(colnames(table),expected))
        {
                stop('The input file is not properly formated:
                -please use the original output excel table from the Connectivity Map (CMap)
                -first-sheet column names shoud be: \n\t\t',paste(expected,collapse=', '))
        }
        rank <- rep(1, dim(table)[1])
        df <- cbind(table[2:dim(table)[2]], rank)
    }

    # pvalues to numeric, na removed afterwards
    df$p <- suppressWarnings(as.numeric(df$p))
    df <- df[-which(df$p > 0.05 | is.na(df$p)), ]

    #cutoff on replicates number
    if (!is.null(n.rep.cutoff)) {
        if (n.rep.cutoff > 1 & n.rep.cutoff > min(df$n)) {
            df <- df[-which(df$n < n.rep.cutoff), ]
        } else if (n.rep.cutoff <= 1) {
            stop("cutoff=NULL already takes all the batches of 1 replicate
            please set a higher cutoff")
        }
    }

    #enrichment type
    if (enrichment == "positive") {
        df <- df[-which(df$enrichment < 0), ]
    } else {
        df <- df[-which(df$enrichment > 0), ]
    }

    #cutoff on enrichment value
    if (!is.null(abs.enrich.cutoff)) {
        if (abs.enrich.cutoff >= 1 | abs.enrich.cutoff <= 0) {
            stop("abs.enrich.cutoff must be set in the interval 0,1")
        } else if (abs.enrich.cutoff > min(abs(df$enrichment))) {
            df <- df[-which(abs(df$enrichment) < abs.enrich.cutoff), ]
        }
    }

    # retreiving the number of cell lines from remaining molecules
    if (plot == "cell.lines") {
        n.cell.lines <- length(unique(df$cell.line))

        i <- 2
        while (i <= n.cell.lines) {
            df$rank[which(df$cell.line == unique(df$cell.line)[i])] <- i
            i <- i + 1
        }
        #max cell line number - colors
        colors <- c("#00AFBB", "#E7B800", "#FC4E07","#81EE81","#DA70D6")
    }

    #jittering parameter
    if (jittering == TRUE) {
        df2 <- df %>% mutate(enrichment2 = jitter(df$enrichment, amount = 0.05))
    } else {
        df2 <- df
        df2$enrichment2 <- df2$enrichment
        df2$rank2 <- df2$rank
    }

    # renaming data
    names(df2)[which(names(df2) == "p")] <- "-log2(pValue)"
    df2$`-log2(pValue)` <- as.numeric(df2$`-log2(pValue)`)
    df2$`-log2(pValue)` <- -log2(df2$`-log2(pValue)` + 0.01)

    #alpha parameter
    alpha <- abs(rescale(df2$n, to = c(0, 1)))

    # assembling data
    if (plot == "cell.lines") {
        df3 <- cbind(df2, alpha, colors[df$rank])
        names(df3)[which(names(df3) == "colors[df$rank]")] <- "cols"
    } else {
        df3 <- cbind(df2, alpha)
    }

    # x final position
    df3$rank2 <- df3$rank + rescale(as.numeric(df3$specificity), to = c(-0.5, 0.5))

    #transparency settings
    breaks.alpha <- c(seq(min(df3$n), max(df3$n), (max(df3$n) / 4)), max(df3$n))
    if (sum(duplicated(breaks.alpha)) > 0) {
        breaks.alpha <- unique(breaks.alpha)
    }
    p <- df3$`-log2(pValue)`
    breaks.p <- 1 / (2^seq(-log2(0.05), max(p), 0.5))
    max.es <- round(max(abs(df3$enrichment)), 2)
    min.es <- round(min(abs(df3$enrichment)), 2)

    #plot settings
    theme_set(
        theme_bw() +
            theme(
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), legend.position = "bottom",
                legend.box = "vertical"
            ) + theme_void()
    )

    #ggplot non global variables - cran note
    rank2 <- enrichment2 <- cols <- drug <- NULL

    #plots
    if (plot == "cell.lines") {
        bubble_pl <- ggplot(df3, aes(x = rank2, y = abs(enrichment2))) +
        geom_point(aes(color = df3$cols, size = df3$`-log2(pValue)`, alpha = alpha)) +
        scale_color_manual(values = colors[seq_len(n.cell.lines)], name = "Cell lines", labels = unique(df3$cell.line)) +
        geom_vline(xintercept = seq_len(n.cell.lines), colour = colors[seq_len(n.cell.lines)], linetype = "dashed", alpha = 0.5) +
        scale_size(range = c(0.5, 10), name = "pValue", labels = round(breaks.p, 2), limits = c(-1, 10)) + # Adjust the range of points size
        scale_alpha("Number of replicats", labels = breaks.alpha) +
        scale_x_continuous(breaks = seq_len(3)) +
        geom_text(aes(label = drug), hjust = 0.5, vjust = 2, alpha = 0.5) +
        annotate("label", x = min(df3$rank2) - 0.1, y = max(abs(df3$enrichment2)) + 0.1, label = paste("ES max", max.es, sep = " "), fill = "grey95") +
        annotate("label", x = min(df3$rank2) - 0.1, y = min(abs(df3$enrichment2)) - 0.1, label = paste("ES min", min.es, sep = " "), fill = "grey95")
    } else {
        bubble_pl <- ggplot(df3, aes(x = rank2, y = abs(enrichment2))) +
        geom_point(aes(color = "#00AFBB", size = df3$`-log2(pValue)`, alpha = alpha)) +
        geom_vline(xintercept = 1, colour = "#FC4E07", linetype = "dashed", alpha = 0.5) +
        geom_text(aes(label = df3$`cmap name`), hjust = 0.5, vjust = 2, alpha = 0.5) +
        scale_size(range = c(0.5, 10), name = "pValue", labels = round(breaks.p, 2), limits = c(-1, 10)) + # Adjust the range of points size
        scale_alpha("Number of replicats", labels = breaks.alpha) +
        annotate("label", x = min(df3$rank2) - 0.1, y = max(abs(df3$enrichment2)) + 0.1, label = paste("ES max", max.es, sep = " "), fill = "grey95") +
        annotate("label", x = min(df3$rank2) - 0.1, y = min(abs(df3$enrichment2)) - 0.1, label = paste("ES min", min.es, sep = " "), fill = "grey95") +
        guides(color = FALSE)
    }

    #output
    if (!(is.null(output_path))) {
        if  (str_sub(output_path, start= -1)!="/")
            {output_path=paste(output_path,"/",sep="")}
        ggsave(paste(output_path,"bubble_cmap.png",sep=""))
        save.image(file = paste(output_path,"bubble.cmap.Rdata",sep=""))
    }

    if (return.gg.table == TRUE) {
        return(df3)
    } else {
        print(bubble_pl)
    }
}
