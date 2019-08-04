CC = gcc
CFLAGS = -c -Wall -march=armv7-a -mfpu=neon -mfloat-abi=hard -mcpu=cortex-a8 -mtune=cortex-a8 -ftree-vectorize -mvectorize-with-neon-quad -O3 -g  -pg -ffast-math 

sg: CoEst_offline.o
	$(CC) CoEst_offline.o -lrt -lm -g -static -o $@  -pg

run:
	$(CC) CoEst_offline.c -o CoEst_offline -lm

clean:
	rm -f CoEst_offline *.o CoEst_offline *.s out_est_SOC.txt sg *.o sg *.s