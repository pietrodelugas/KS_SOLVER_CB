
#include <stdio.h>
#include <sys/resource.h>

int subcheckmem() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    printf("Max RSS: %ld %s\n", ru.ru_maxrss, "KB");
    return 0;
}
