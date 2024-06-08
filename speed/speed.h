#ifndef SPEED_H
#define SPEED_H

#define DoNotOptimize(value) asm volatile("" : : "r,m"(value) : "memory")

#endif  // SPEED_H
