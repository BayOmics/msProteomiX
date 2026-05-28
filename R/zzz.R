# ==============================================================================
# msProteomiX — 包加载钩子
# ==============================================================================

.onLoad <- function(libname, pkgname) {
  # .msProteomiX_env 已在 utils.R 中定义，此处无需重新创建
  # 仅在需要时清理旧数据
  if (exists(".msProteomiX_env", envir = parent.env(environment()))) {
    env <- get(".msProteomiX_env", envir = parent.env(environment()))
    rm(list = ls(env), envir = env)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "=== msProteomiX v", utils::packageVersion(pkgname), " ===\n",
    "Multi-Engine Mass Spectrometry Proteomics Data Processing\n",
    "Supported engines: FragPipe | Spectronaut | MaxQuant | PD | DIA-NN\n",
    "Type ?read_ms_data to get started."
  )
}
