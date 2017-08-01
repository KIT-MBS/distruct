
/*************************************
 *
 *  Filename : sig.h
 *
 *  Projectname : 
 *
 *  Author : Oskar Taubert
 *
 *  Creation Date : Mon 31 Jul 2017 11:05:33 CEST
 *
 *  Last Modified : Mon 31 Jul 2017 11:25:14 CEST
 *
 * *************************************/

#ifndef SIG_H
#define SIG_H

#include <execinfo.h>
#include <signal.h>
//#include <cstdio>
#include <stdio.h>
#include <stdlib.h>

static void handle_signals(int signal)
{
    switch(signal)
    {
        case SIGINT:
            {
                fputs("received SIGINT keyboard interrupt.\n", stderr);
            }
            break;
        // TODO cases SIGABRT, SIGFPE, SIGILL, SIGSEGV
        case SIGTERM:
            {
            }
        default:
            {
                fputs("received SIGTERM termination signal.\n", stderr);
            }
            break;
    }
    exit(signal);
}

#endif /* SIG_H */
