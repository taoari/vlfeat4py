SRC=../src
INC=../include

PROG=gmm
g++ -o test_vl_$PROG $SRC/vl_$PROG.cpp $SRC/vl_dsift.cpp test_vl_$PROG.cpp -I $INC -larmadillo -lvl
./test_vl_$PROG
# valgrind --leak-check=summary ./test_vl_$PROG
