# ============================================================
# enrichment_bubbleplot.R —— 通用 GO 富集气泡图（可复用）
#
# 数据不写死在脚本里：每次从输入文件读取，换数据只换输入文件。
#
# 两个输入：
#   1) 绘图数据（CSV/TSV）：手动筛选后要画的 term，固定 5 列：
#        term_id        → GO ID（仅存档，不绘制）
#        description    → Y 轴描述
#        count          → 气泡大小
#        strength       → X 轴（脚本内做 10^strength）
#        fdr            → 颜色（每个 term 自己的 -log10(fdr)）
#   2) 原始富集结果（CSV/TSV）：只用来抓色标范围——
#        色标上下限 = 原始结果全部 term 的 -log10(FDR) 最小/最大，
#        保证换数据 / 换筛选后色标一致，不随当前画出的 term 缩放。
#
# 筛选：FDR 自然断点——仅当绘图数据 term 数 > 15 时启用；
#       相邻 -log10(FDR) 落差最大处断开，取上方更显著的一簇。
# 输出：175mm 宽 × 自适应高、900dpi、LZW 的 TIFF
#
# 用法：
#   Rscript enrichment_bubbleplot.R 绘图数据.csv 原始结果.tsv [输出.tiff]
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# 字体：Windows 用 Arial；其它系统退回默认 sans（避免因缺 Arial 直接报错）
FONT <- if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  "Arial"
} else {
  "sans"
}

# ---- 命令行参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("用法：Rscript enrichment_bubbleplot.R <绘图数据.csv|tsv> <原始结果.csv|tsv> [输出.tiff]")
in_plot <- args[1]
in_raw  <- args[2]
if (!file.exists(in_plot)) stop("找不到绘图数据文件：", in_plot)
if (!file.exists(in_raw))  stop("找不到原始结果文件：", in_raw)
out_path <- if (length(args) >= 3 && nzchar(args[3])) args[3] else
  file.path(dirname(in_plot),
            paste0(sub("\\.[^.]+$", "", basename(in_plot)), "_bubbleplot_900dpi.tiff"))

# ---- 读表格（TSV / CSV 自动识别）----
read_table <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  if (length(first) == 0 || !nzchar(trimws(first)))
    stop("输入文件是空的，没有读到任何内容：", path)
  raw <- if (grepl("\t", first)) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  }
  if (nrow(raw) == 0)
    stop("这个文件只有表头、没有任何数据行：", path)
  raw
}

# ---- 1) 绘图数据（固定 5 列）----
raw_plot <- read_table(in_plot)
get_col <- function(name) {
  i <- match(name, names(raw_plot))
  if (is.na(i)) stop("绘图数据缺少必需列「", name, "」，请检查表头。当前文件里的列是：", paste(names(raw_plot), collapse = " / "))
  raw_plot[[i]]
}
d_full <- data.frame(
  term_id     = as.character(get_col("term_id")),
  description = as.character(get_col("description")),
  count       = as.numeric(get_col("count")),
  strength    = as.numeric(get_col("strength")),
  fdr         = as.numeric(get_col("fdr")),
  stringsAsFactors = FALSE
)
d_full <- d_full[!is.na(d_full$count) & !is.na(d_full$strength) &
                 !is.na(d_full$fdr) & d_full$fdr > 0 & nzchar(d_full$description), ]
if (nrow(d_full) == 0)
  stop("绘图数据清洗后没有留下有效行（可能是必需列的值缺失或格式不对）。请检查输入数据。")
if (nrow(d_full) == 1)
  stop("只有 1 个 term，画不出有意义的气泡图。请至少提供 2 个 term。")
if (all(d_full$fdr > 0.05))
  stop("所有 term 的 FDR 都大于 0.05，没有统计学显著的富集结果。建议重新进行富集分析（或调整参数/背景基因集）后再出图。")
d_full <- d_full[order(d_full$fdr), ]   # FDR 升序（最显著在顶部）

# ---- 2) 原始富集结果：只抓色标范围 ----
find_fdr_col <- function(nms) {
  n <- tolower(trimws(nms))
  for (p in c("false discovery rate", "fdr", "q-value", "qvalue",
              "adjusted p", "corrected p", "p.adjust", "padj")) {
    i <- which(n == p | grepl(p, n, fixed = TRUE))
    if (length(i)) return(i[1])
  }
  NA
}
raw_raw <- read_table(in_raw)
i_fdr <- find_fdr_col(names(raw_raw))
if (is.na(i_fdr))
  stop("原始结果里找不到 FDR 列（false discovery rate / fdr / q-value 等）。当前列是：", paste(names(raw_raw), collapse = " / "))
raw_fdr <- suppressWarnings(as.numeric(raw_raw[[i_fdr]]))
raw_fdr <- raw_fdr[!is.na(raw_fdr) & raw_fdr > 0]
if (length(raw_fdr) == 0)
  stop("原始结果的 FDR 列没有有效数值（应为 0~1 之间、大于 0 的数）。")
col_lo <- min(-log10(raw_fdr))
col_hi <- max(-log10(raw_fdr))

# ---- 筛选：FDR 自然断点（仅当绘图数据 > 15 个 term 时启用）----
NATURAL_BREAK_MIN <- 15
n_all <- nrow(d_full)
break_at <- NA
d <- d_full
if (n_all > NATURAL_BREAK_MIN) {
  all_nl <- -log10(d_full$fdr)          # 已按 FDR 升序
  gaps <- diff(all_nl)                  # 相邻落差（均为负，越负落差越大）
  k <- which.min(gaps)
  if (gaps[k] < 0) {                    # 存在真实落差才断开
    d <- d_full[seq_len(k), ]
    break_at <- k
  }
}

# ---- 绘图变量 ----
d$negLog10FDR <- -log10(d$fdr)
d$RichFactor  <- 10^d$strength

# ---- 自适应：X 轴刻度（最粗步长 100→50→10，刻度数 3~6，覆盖最大值）----
max_fold <- max(d$RichFactor)
xbreaks <- NULL
for (st in c(100, 50, 10)) {
  b <- seq(0, ceiling(max_fold / st) * st, by = st)
  if (length(b) >= 3 && length(b) <= 6) { xbreaks <- b; break }
}
if (is.null(xbreaks)) xbreaks <- seq(0, ceiling(max_fold / 100) * 100, by = 100)
xmax <- max(xbreaks)

# ---- 自适应：尺寸图例 count 值（min~max，直径线性映射 1.0~3.0mm）----
cnt_lo <- min(d$count); cnt_hi <- max(d$count)
sz_counts <- seq(cnt_lo, cnt_hi, by = 1)
sz_diam   <- 1.0 + (sz_counts - cnt_lo) / max(1, cnt_hi - cnt_lo) * 2.0
sz_radius <- sz_diam / 2

# ---- 自适应：画布高度 + 图例垂直居中 ----
n <- nrow(d)
height_mm <- min(110, max(70, 40 + n * 8))
shift_y   <- height_mm / 2 - 34.4

# ---- 图例几何坐标（水平固定）----
leg_x_class  <- 102.35;  leg_y_class  <- 30.4 + shift_y
leg_x_bubble <- 103.65;  leg_y_bubble <- 33.6 + shift_y
leg_x_go     <- 106.55;  leg_y_go     <- 33.6 + shift_y
bubble_r     <- 0.5
sz_bubble_x  <- 103.55
sz_y         <- 43.0 + shift_y + (0:(length(sz_counts) - 1)) * 3.06
sz_num_x     <- 105.85
sz_label_x   <- 102.25
sz_label_y   <- 39.7 + shift_y

# 描述超长换行
wrap_label <- function(x, width = 70)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
d$description <- wrap_label(d$description)
d$description <- factor(d$description, levels = rev(d$description))

# ---- 主图 ----
p <- ggplot(d, aes(x = RichFactor, y = description)) +
  geom_point(aes(size = count, colour = negLog10FDR)) +
  scale_colour_gradientn(
    colours = c("#00FF00", "#94D900", "#BEB700", "#F06100", "#FD2100", "#FF0000"),
    name = expression(-log[10](pvalue)),
    limits = c(col_lo, col_hi),
    guide = guide_colourbar(barwidth = unit(3, "mm"), barheight = unit(15.4, "mm"), order = 1)
  ) +
  scale_size_continuous(range = c(1.0, 3.0), guide = "none") +
  scale_x_continuous(breaks = xbreaks, limits = c(0, xmax)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8, base_family = FONT) +
  theme(
    text             = element_text(family = FONT),
    axis.text        = element_text(colour = "black", family = FONT),
    axis.title       = element_text(family = FONT),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    legend.title     = element_text(size = 6.4, family = FONT),
    legend.text      = element_text(size = 6.4, family = FONT),
    legend.justification = "top",
    legend.margin    = margin(t = 4.2 + shift_y, r = 0, b = 0, l = 0, unit = "mm"),
    plot.margin      = margin(t = 2, r = 58, b = 2, l = 2, unit = "mm")
  )

# ---- 输出 ----
tiff(out_path, width = 175, height = height_mm, units = "mm", res = 900, compression = "lzw")
gt <- ggplotGrob(p)
grid.newpage()
grid.draw(gt)

# 「class + 一个气泡 + GO」单元
grid.text("class", x = unit(leg_x_class/175, "npc"), y = unit((height_mm - leg_y_class)/height_mm, "npc"),
          just = c("left", "center"), gp = gpar(fontsize = 6.4, fontfamily = FONT, col = "black"))
grid.circle(x = unit(leg_x_bubble/175, "npc"), y = unit((height_mm - leg_y_bubble)/height_mm, "npc"),
            r = unit(bubble_r, "mm"), gp = gpar(fill = "grey30", col = NA))
grid.text("GO", x = unit(leg_x_go/175, "npc"), y = unit((height_mm - leg_y_go)/height_mm, "npc"),
          just = c("left", "center"), gp = gpar(fontsize = 6.4, fontfamily = FONT, col = "black"))

# 「count + N 个气泡 + 数字」尺寸图例
grid.text("count", x = unit(sz_label_x/175, "npc"), y = unit((height_mm - sz_label_y)/height_mm, "npc"),
          just = c("left", "center"), gp = gpar(fontsize = 6.4, fontfamily = FONT, col = "black"))
for (i in seq_along(sz_counts)) {
  grid.circle(x = unit(sz_bubble_x/175, "npc"), y = unit((height_mm - sz_y[i])/height_mm, "npc"),
              r = unit(sz_radius[i], "mm"), gp = gpar(fill = "black", col = NA))
  grid.text(as.character(sz_counts[i]), x = unit(sz_num_x/175, "npc"), y = unit((height_mm - sz_y[i])/height_mm, "npc"),
            just = c("left", "center"), gp = gpar(fontsize = 6.4, fontfamily = FONT, col = "black"))
}
dev.off()

# ---- 摘要 ----
cat("已输出:", out_path, "\n")
cat("绘图数据 term 数:", n_all, "→ 自然断点后:", n, "\n")
if (!is.na(break_at))
  cat("断点位置: FDR", format(d_full$fdr[break_at], scientific = TRUE), "→",
      format(d_full$fdr[break_at + 1], scientific = TRUE), "\n")
cat("色标 -log10(FDR)（来自原始结果", length(raw_fdr), "个 term）:", round(col_lo, 3), "→", round(col_hi, 3), "\n")
cat("X 刻度:", paste(xbreaks, collapse = " / "), "\n")
cat("富集倍数范围:", round(min(d$RichFactor), 2), "→", round(max(d$RichFactor), 2), "\n")
