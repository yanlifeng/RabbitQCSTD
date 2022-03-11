#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "stats.h"
#include "writer.h"
#include "options.h"
#include "filterresult.h"
#include "TGSStats.h"

using namespace std;

class ThreadConfig {
public:
    ThreadConfig(Options *opt, int threadId, bool paired = false);

    ~ThreadConfig();

    inline Stats *getPreStats1() { return mPreStats1; }

    inline Stats *getPostStats1() { return mPostStats1; }

    inline Stats *getPreStats2() { return mPreStats2; }

    inline Stats *getPostStats2() { return mPostStats2; }

    inline TGSStats *getTGSStats() { return mTGSStats; }

    inline Writer *getWriter1() { return mWriter1; }

    inline Writer *getWriter2() { return mWriter2; }

    inline FilterResult *getFilterResult() { return mFilterResult; }

    void initWriter(string filename1);

    void initWriter(string filename1, string filename2);

    void initWriter(ofstream *stream);

    void initWriter(ofstream *stream1, ofstream *stream2);

    void initWriter(gzFile gzfile);

    void initWriter(gzFile gzfile1, gzFile gzfile2);

    void addFilterResult(int result);

    int getThreadId() { return mThreadId; }

    // for splitting output
    // increase mCurrentSplitReads by readNum, and check it with options->split.size;
    void markProcessed(long readNum);

    void initWriterForSplit();

    bool canBeStopped();

    void cleanup();

private:
    void deleteWriter();

    void writeEmptyFilesForSplitting();

public:
    double cost;
    double cost1;
    double cost2;
    double cost3;
    double cost4;
    double cost5;
    double cost6;
    double cost7;
    double cost8;
    double cost9;
    double cost10;
    double cost11;
    double cost12;
    double cost13;
    double cost14;
    double costFormat;
    int totCnt;
    Stats *mPreStats1;
    Stats *mPostStats1;
    Stats *mPreStats2;
    Stats *mPostStats2;
    Writer *mWriter1;
    Writer *mWriter2;
    Options *mOptions;
    FilterResult *mFilterResult;

    // for spliting output
    int mThreadId;
    int mWorkingSplit;
    long mCurrentSplitReads;
    bool mCanBeStopped;
    //---for tgs---//
    TGSStats *mTGSStats;
};

#endif
