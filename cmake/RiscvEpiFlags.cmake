set(QE_EPI_COMPILE_OPTIONS -O3 -mcpu=atrevido-vec -ffast-math -mepi -mllvm -combiner-store-merging=0 -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -mllvm -vectorizer-use-vp-strided-load-store -mllvm -disable-loop-idiom-memcpy -mllvm -disable-loop-idiom-memset -Rpass-missed=loop-vectorize -Xflang -target-feature -Xflang +does-not-implement-vszext -Xflang -target-feature -Xflang +does-not-implement-tu -mllvm -riscv-uleb128-reloc=0 -fno-slp-vectorize)
message(STATUS "set epi compile options to ${QE_EPI_COMPILE_OPTIONS}")

