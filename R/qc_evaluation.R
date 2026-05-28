# ==============================================================================
# msProteomiX — 综合评估 QC 模块
# ==============================================================================

#' 计算等电点 pI (Bjellqvist scale)
#'
#' @param seq 氨基酸序列字符串
#' @return pI 数值
#' @export
#' @examples
#' calc_pI("ACDEFGHIKLMNPQRSTVWY")
calc_pI <- function(seq) {
  if (is.na(seq) || seq == "") return(NA_real_)
  seq <- toupper(gsub("[^A-Z]", "", seq))
  aa <- unlist(strsplit(seq, ""))

  # pKa values (Bjellqvist scale)
  pKa <- c(N_term = 9.6, K = 10.0, R = 12.0, H = 5.98,
           D = 4.05, E = 4.45, C = 9.0, Y = 10.0, C_term = 2.36)

  cnts <- table(aa)
  .get <- function(x) { v <- as.integer(cnts[x]); ifelse(is.na(v), 0L, v) }

  count_K <- .get("K"); count_R <- .get("R"); count_H <- .get("H")
  count_D <- .get("D"); count_E <- .get("E"); count_C <- .get("C")
  count_Y <- .get("Y")

  charge <- function(pH) {
    pos <- (10^pKa["N_term"] / (10^pKa["N_term"] + 10^pH)) +
      (count_K * 10^pKa["K"] / (10^pKa["K"] + 10^pH)) +
      (count_R * 10^pKa["R"] / (10^pKa["R"] + 10^pH)) +
      (count_H * 10^pKa["H"] / (10^pKa["H"] + 10^pH))

    neg <- (10^pH / (10^pKa["C_term"] + 10^pH)) +
      (count_D * 10^pH / (10^pKa["D"] + 10^pH)) +
      (count_E * 10^pH / (10^pKa["E"] + 10^pH)) +
      (count_C * 10^pH / (10^pKa["C"] + 10^pH)) +
      (count_Y * 10^pH / (10^pKa["Y"] + 10^pH))
    pos - neg
  }

  min_pH <- 0; max_pH <- 14
  for (i in 1:15) {
    mid <- (min_pH + max_pH) / 2
    if (charge(mid) > 0) min_pH <- mid else max_pH <- mid
  }
  (min_pH + max_pH) / 2
}


#' 计算 GRAVY 疏水性指数
#'
#' @param seqs 氨基酸序列字符串向量
#' @return GRAVY 数值向量
#' @export
calc_gravy <- function(seqs) {
  hydropathy <- c(A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5,
                  Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I = 4.5,
                  L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6,
                  S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2)
  sapply(seqs, function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    s_clean <- gsub("[^A-Z]", "", toupper(s))
    aa <- unlist(strsplit(s_clean, ""))
    scores <- hydropathy[aa]
    if (length(scores) == 0 || all(is.na(scores))) return(NA_real_)
    mean(scores, na.rm = TRUE)
  }, USE.NAMES = FALSE)
}


#' 计算 Missed Cleavage 数量
#'
#' @param seqs 肽段序列字符串向量
#' @return 整数向量
#' @export
calc_missed_cleavage <- function(seqs) {
  sapply(seqs, function(s) {
    if (is.na(s) || nchar(s) <= 1) return(0L)
    int_seq <- substr(s, 1, nchar(s) - 1)
    stringr::str_count(int_seq, "[KR](?!P)")
  }, USE.NAMES = FALSE)
}
