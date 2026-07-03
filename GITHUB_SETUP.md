# scONCO — 獨立 package repo (根目錄就是 package)

這個資料夾的**根目錄就是 R package**,所以 `install_github` **不需要 subdir**,也符合投稿與一般 R package 慣例。

已內建 (self-contained):
- `R/` — 對外函數 + roxygen 註解
- `inst/legacy/R/` — 13 個核心演算法模組
- `inst/extdata/` — marker DB (human/mouse + 各版本)

## 一、改一個地方
編輯 `DESCRIPTION`,把 `YOUR-GITHUB-ACCOUNT` 換成你的 GitHub 帳號 (兩行):
```
URL: https://github.com/<你的帳號>/scONCO
BugReports: https://github.com/<你的帳號>/scONCO/issues
```

## 二、產生說明文件 (讓 man/ 不再空)
```r
install.packages(c("devtools","roxygen2"))
devtools::document(".")   # 讀 R/ 的 roxygen 註解 → 生成 man/*.Rd + 更新 NAMESPACE
```

## 三、本機自我檢查 (投稿前)
```r
devtools::check(".")      # R CMD check;把 error/warning 貼給我,我幫你修
```

## 四、在 GitHub 建空 repo → 推上去
到 github.com 建一個名為 `scONCO` 的空 repo (不勾 README/.gitignore),然後在這個資料夾:
```bash
cd E:/Charlene/Bioinformatics_Tool_Development/scONCO_tool
git init
git add .
git commit -m "scONCO v1.1: pan-cancer scRNA-seq annotation package"
git branch -M main
git remote add origin https://github.com/<你的帳號>/scONCO.git
git push -u origin main
```
> push 要密碼時,用 GitHub Personal Access Token (Settings → Developer settings → Tokens),不是登入密碼。

## 五、別人怎麼安裝 (乾淨,不用 subdir)
```r
remotes::install_github("<你的帳號>/scONCO")
library(scONCO)
markers_df <- load_cancer_marker_db()
# seu <- run_scONCO(seu, markers_df, apply_confidence = TRUE)
```

## 與大 repo 的關係
- `scONCO/`  (原本的) = **研究專案**:benchmarks、manuscript、database、圖表、所有分析
- `scONCO_tool/` (這個) = **發行用 package**:只有 package 本體,乾淨、可 install_github

兩者分開是 method paper 常見做法:大 repo 給 reviewer 看完整分析、拿 Zenodo DOI;package repo 給使用者一鍵安裝。DB 與核心已複製進 package,兩邊獨立不互相依賴。
