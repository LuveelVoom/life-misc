gcc transcensus.c -lm -o tc128.exe -std=c99 -O3 -Wall -Wextra -Wno-unused-function -Werror -fno-tree-loop-distribute-patterns -march=native -D __NO_AVX2
gcc transcensus.c -lm -o tc256.exe -std=c99 -O3 -Wall -Wextra -Wno-unused-function -Werror -fno-tree-loop-distribute-patterns -march=native
