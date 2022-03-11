#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"


SingleEndProcessor::SingleEndProcessor(Options *opt) {
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream = NULL;
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter = NULL;

    mDuplicate = NULL;
    if (mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    //dsrc faster reader
    fastqPool = new dsrc::fq::FastqDataPool(128, SwapBufferSize);
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if (mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    //delete dsrc mem pool
    delete fastqPool;
}

void SingleEndProcessor::initOutput() {
    if (mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if (mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig *config) {
    if (mOptions->out1.empty())
        return;

    if (mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool SingleEndProcessor::process() {
    if (!mOptions->split.enabled)
        initOutput();

    initPackRepository();
    std::thread
            producer(std::bind(&SingleEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig **configs = new ThreadConfig *[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
    }

    std::thread **threads = new thread *[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread *leftWriterThread = NULL;
    if (mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));

    printf("all thread created\n");
    producer.join();
    printf("producer join\n");
    for (int t = 0; t < mOptions->thread; t++) {
        threads[t]->join();
    }
    printf("consumers join\n");

    if (!mOptions->split.enabled) {
        if (leftWriterThread)
            leftWriterThread->join();
        printf("writer join\n");

    }

    if (mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats *> preStats;
    vector<Stats *> postStats;
    vector<FilterResult *> filterResults;
    for (int t = 0; t < mOptions->thread; t++) {
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats *finalPreStats = Stats::merge(preStats);
    Stats *finalPostStats = Stats::merge(postStats);
    FilterResult *finalFilterResult = FilterResult::merge(filterResults);


    auto OverRepSeq1 = finalPreStats->mOverRepSeq;
//    cout << "=============OverRepSeq1=============" << endl;
//    for (auto it:OverRepSeq1) {
//        cout << it.first << " " << it.second << endl;
//    }
//    cout << "=====================================" << endl;

    auto OverRepSeq2 = finalPostStats->mOverRepSeq;
//    cout << "=============OverRepSeq2=============" << endl;
//    for (auto it:OverRepSeq2) {
//        cout << it.first << " " << it.second << endl;
//    }
//    cout << "=====================================" << endl;


    ofstream ofs;
    ofs.open("ORP2.log", ifstream::out);
    for (auto it:OverRepSeq1) {
        ofs << it.first << " " << it.second << "\n";
    }
    ofs.close();

    ofs.open("ORP3.log", ifstream::out);
    for (auto it:OverRepSeq2) {
        ofs << it.first << " " << it.second << "\n";
    }
    ofs.close();

//
//    printf("after size %d\n", finalPostStats->getCycles());
//    auto p1 = finalPostStats->mCycleBaseContents;
//    auto p2 = finalPostStats->mCycleBaseQual;
//    auto p3 = finalPostStats->mQualityCurves["mean"];
//    for (int i = 0; i < finalPostStats->getCycles(); i++) {
//
//        if (i > 900) {
//            printf("i ========== %d\n", i);
//
//            int64_t sum_qul = 0;
//            int64_t sum_cnt = 0;
//            for (int j = 0; j < 8; j++) {
//                sum_qul += p2[j][i];
//                sum_cnt += p1[j][i];
//            }
//            if (sum_cnt == 0) {
//                printf("sum_cnt is 0\n");
//                continue;
//            }
//            double res = 1.0 * sum_qul / sum_cnt;
//            printf("sum_qul %lld , sum_cnt %lld\n", sum_qul, sum_cnt);
//            for (int j = 0; j < 8; j++)printf("%ld ", p1[j][i]);
//            printf("\n");
//            for (int j = 0; j < 8; j++)printf("%ld ", p2[j][i]);
//            printf("\n");
//            printf("tmp_double %lf\n", res);
//            printf("rrrr %lf\n", p3[i]);
//        }
//    }

    // read filter results to the first thread's
//    for (int t = 1; t < mOptions->thread; t++) {
//        preStats.push_back(configs[t]->getPreStats1());
//        postStats.push_back(configs[t]->getPostStats1());
//    }

//    fstream out;
//
//    string outFileName1 = "preStatsKmer";
//    string outFileName2 = "postStatsKmer";
//    string outFileName3 = "mDuplicateCount";
//    string outFileName4 = "mDuplicateGC";
//    string outFileName5 = "mDuplicateDups";
//    string outFileName6 = "preStatsTotalBase";
//    string outFileName7 = "preStatsTotalQual";
//    string outFileName8 = "postStatsTotalBase";
//    string outFileName9 = "postStatsTotalQual";
//
//    out.open(outFileName1.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPreStats->mKmer), finalPreStats->mKmerBufLen * sizeof(long));
//    out.close();
//
//    out.open(outFileName2.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPostStats->mKmer), finalPostStats->mKmerBufLen * sizeof(long));
//    out.close();
//
//    out.open(outFileName3.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(mDuplicate->mCounts), mDuplicate->mKeyLenInBit * sizeof(uint16));
//    out.close();
//
//    out.open(outFileName4.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(mDuplicate->mGC), mDuplicate->mKeyLenInBit * sizeof(uint8));
//    out.close();
//
//    out.open(outFileName5.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(mDuplicate->mDups), mDuplicate->mKeyLenInBit * sizeof(uint64));
//    out.close();
//
//    out.open(outFileName6.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPreStats->mCycleTotalBase), finalPreStats->mBufLen * sizeof(long));
//    out.close();
//
//    out.open(outFileName7.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPreStats->mCycleTotalQual), finalPreStats->mBufLen * sizeof(long));
//    out.close();
//
//    out.open(outFileName8.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPostStats->mCycleTotalBase), finalPostStats->mBufLen * sizeof(long));
//    out.close();
//
//    out.open(outFileName9.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(finalPostStats->mCycleTotalQual), finalPostStats->mBufLen * sizeof(long));
//    out.close();

#ifdef Timer
    double cost = 0;
    double cost1 = 0;
    double cost2 = 0;
    double cost3 = 0;
    double cost4 = 0;
    double cost5 = 0;
    double cost6 = 0;
    double cost7 = 0;
    double cost8 = 0;
    double cost9 = 0;
    double cost10 = 0;
    double cost11 = 0;
    double cost12 = 0;
    double cost13 = 0;
    double cost14 = 0;
    double costFormat = 0;
    int totCnt = 0;

    for (int t = 0; t < mOptions->thread; t++) {
        cost += configs[t]->cost;
        cost1 += configs[t]->cost1;
        cost2 += configs[t]->cost2;
        cost3 += configs[t]->cost3;
        cost4 += configs[t]->cost4;
        cost5 += configs[t]->cost5;
        cost6 += configs[t]->cost6;
        cost7 += configs[t]->cost7;
        cost8 += configs[t]->cost8;
        cost9 += configs[t]->cost9;
        cost10 += configs[t]->cost10;
        cost11 += configs[t]->cost11;
        cost12 += configs[t]->cost12;
        cost13 += configs[t]->cost13;
        cost14 += configs[t]->cost14;
        totCnt += configs[t]->totCnt;
        costFormat += configs[t]->costFormat;

    }


    printf("total getPreStats1()->statRead(or1) ====: %.5f\n", cost1);
    printf("total mDuplicate->statRead(or1) ========: %.5f\n", cost2);
    printf("total mOptions->indexFilter()  =========: %.5f\n", cost3);
    printf("total mUmiProcessor->process(or1) ======: %.5f\n", cost4);
    printf("total mFilter->trimAndCut() ============: %.5f\n", cost5);
    printf("total PolyX::trimPolyG() ===============: %.5f\n", cost4);
    printf("total trimBySequence ===================: %.5f\n", cost7);
    printf("total r1->resize() =====================: %.5f\n", cost8);
    printf("total mFilter->passFilter(r1) ==========: %.5f\n", cost9);
    printf("total addFilterResult(result) ==========: %.5f\n", cost10);
    printf("total outstr += r1->toString() =========: %.5f\n", cost11);
    printf("total getPostStats1()->statRead(r1) ====: %.5f\n", cost12);
    printf("total delete r1 ========================: %.5f\n", cost13);
    printf("total ready output ========================: %.5f\n", cost14);
    printf("total costTotel ========================: %.5f\n",
           cost1 + cost2 + cost3 + cost4 + cost5 + cost6 + cost7 + cost8 + cost9 + cost10 + cost11 + cost12 + cost13);
    printf("total cost =============================: %.5f\n", cost);
    printf("total  =================================: %d\n", totCnt);
    printf("total format =================================: %.5f\n", costFormat);

#endif
    cerr << "Read1 before filtering:" << endl;
    finalPreStats->print();
    cerr << endl;
    cerr << "Read1 after filtering:" << endl;
    finalPostStats->print();

    cerr << endl;
    cerr << "Filtering result:" << endl;
    finalFilterResult->print();

    int *dupHist = NULL;
    double *dupMeanTlen = NULL;
    double *dupMeanGC = NULL;
    double dupRate = 0.0;
    if (mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for (int t = 0; t < mOptions->thread; t++) {
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if (mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if (leftWriterThread)
        delete leftWriterThread;

    if (!mOptions->split.enabled)
        closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack *pack, ThreadConfig *config) {


    //debug for losing last line
    //cerr << pack->data[pack->count - 1]->mName << endl;
    //cerr << pack->data[pack->count - 1]->mQuality << endl;
    //debug
    string outstr;
    int readPassed = 0;
    //------------------my thinking---------------------------------
    /*
    if (mOptions -> thirdgene){
      for(int p = 0; p < pack->data[p]; p++){
        Read* or1 = pack->data[p];
        config -> getTGStats() -> tgsStatRead(or1);
      }
    }
    */
    //1. add TGStats class in project
    //2. add getTGStats function in ThreadConfig file;
    //3. add TGStats menber variable in ThreadConfig class
    //---------------------------------------------------

#ifdef Timer
    config->totCnt += 1;
    double t0, t1;
    t0 = get_wall_time();
#endif

    for (int p = 0; p < pack->count; p++) {

        // original read1
        Read *or1 = pack->data[p];
        // stats the original read before trimming
#ifdef Timer

        t1 = get_wall_time();
#endif
        config->getPreStats1()->statRead(or1);//cost 12/57
#ifdef Timer

        config->cost1 += get_wall_time() - t1;
        // handling the duplication profiling
        t1 = get_wall_time();
#endif
        //TODO maybe mDuplicate is not thread safe
        if (mDuplicate)
            mDuplicate->statRead(or1);//cost 14/57

#ifdef Timer
        config->cost2 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif

        // filter by index
        if (mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
//            printf("mOptions->indexFilter ...");
            delete or1;
            continue;
        }
#ifdef Timer
        config->cost3 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        // umi processing

        //TODO maybe this can be the big hotspot
        if (mOptions->umi.enabled) {
//            printf("mOptions->umi ...");
            mUmiProcessor->process(or1);
        }
#ifdef Timer
        config->cost4 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        // trim in head and tail, and apply quality cut in sliding window
        //not in
        //TODO look what will do in this function and how to open this function
        //TODO maybe this can be the big hotspot
//        printf("_______________SSSS\n");
//        cout << or1->length() << " " << mOptions->trim.front1 << " " << mOptions->trim.tail1 << endl;
//        printf("_______________DDDD\n");

        Read *r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
//        if (r1 != or1) {
//            cout << r1 << " " << or1 << endl;
//            cout << r1->mSeq.mStr << endl;
//            cout << or1->mSeq.mStr << endl;
//        }
//        if (r1 != NULL && r1 != or1) {
//            printf("GG\n");
//        }
#ifdef Timer
        config->cost5 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        if (r1 != NULL) {
            //not in because xx is false
            if (mOptions->polyGTrim.enabled) {
//                printf("polyGTrim ...\n");
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);

            }
            //not in because xx is false
            if (mOptions->polyXTrim.enabled) {
//                printf("polyXTrim ...\n");
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);

            }
        }
#ifdef Timer
        config->cost6 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        //not in because mOptions->adapter.hasSeqR1 is false
        if (r1 != NULL && mOptions->adapter.enabled && mOptions->adapter.hasSeqR1) {
//            printf("AdapterTrimmer::trimBySequence ...\n");
            AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence);
        }

#ifdef Timer
        config->cost7 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        if (r1 != NULL) {
            // not in because mOptions->trim.maxLen1 is 0
            if (mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length()) {
//                printf("r1->resize ...\n");
                r1->resize(mOptions->trim.maxLen1);
            }
        }
#ifdef Timer
        config->cost8 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        //in ! O(len)
        int result = mFilter->passFilter(r1);
#ifdef Timer
        config->cost9 += get_wall_time() - t1;
        t1 = get_wall_time();
#endif
        config->addFilterResult(result);
#ifdef Timer
        config->cost10 += get_wall_time() - t1;
#endif
        //in !
//        if (result != PASS_FILTER) {
//            cout << r1->toString() << endl;
//        }
        if (r1 != NULL && result == PASS_FILTER) {
#ifdef Timer
            t1 = get_wall_time();
#endif
            outstr += r1->toString();
#ifdef Timer
            config->cost11 += get_wall_time() - t1;
            // stats the read after filtering
            t1 = get_wall_time();
#endif
            config->getPostStats1()->statRead(r1);//cost 23/57
#ifdef Timer
            config->cost12 += get_wall_time() - t1;
#endif
            readPassed++;
        }


//        if (r1 != NULL) {
//            assert(r1->mName == or1->mName);
//            assert(r1->mSeq.mStr == or1->mSeq.mStr);
//            assert(r1->mStrand == or1->mStrand);
//            assert(r1->mQuality == or1->mQuality);
//        }
#ifdef Timer
        t1 = get_wall_time();
#endif
        delete or1;
        // if no trimming applied, r1 should be identical to or1
        //TODO isn't r1 should be identical to or1 in anytime ?
//        if (r1 != or1 && r1 != NULL)
//            delete r1;
#ifdef Timer
        config->cost13 += get_wall_time() - t1;
#endif
    }
#ifdef Timer
    config->cost += get_wall_time() - t0;
    t0 = get_wall_time();
#endif
    // if splitting output, then no lock is need since different threads write different files
    if (!mOptions->split.enabled)
        mOutputMtx.lock();
    if (mOptions->outputToSTDOUT) {
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    } else if (mOptions->split.enabled) {
        // split output by each worker thread
        if (!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr);
    } else {
        if (mLeftWriter) {
            char *ldata = new char[outstr.size()];
            memcpy(ldata, outstr.c_str(), outstr.size());
            mLeftWriter->input(ldata, outstr.size());
        }
    }

    if (!mOptions->split.enabled)
        mOutputMtx.unlock();

    if (mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    //delete pack->data;
    std::vector<Read *>().swap(pack->data);
    delete pack;
#ifdef Timer
    config->cost14 += get_wall_time() - t0;
#endif
    return true;
}


void SingleEndProcessor::initPackRepository() {
    //mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    mRepo.packBuffer = new dsrc::fq::FastqDataChunk *[PACK_NUM_LIMIT];
    //memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    memset(mRepo.packBuffer, 0, sizeof(dsrc::fq::FastqDataChunk *) * PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    //mRepo.readCounter = 0;

}

void SingleEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(dsrc::fq::FastqDataChunk *pack) {
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
    //cerr << "chunk in producePack" << endl;
    //cerr << (char *)pack->data.Pointer() << endl;
    mRepo.writePos++;

    /*if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;*/

    //mRepo.repoNotEmpty.notify_all();
    //lock.unlock();
}

void SingleEndProcessor::consumePack(ThreadConfig *config) {
    dsrc::fq::FastqDataChunk *chunk;
    ReadPack *data = new ReadPack;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    mInputMtx.lock();
    while (mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if (mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    //data = mRepo.packBuffer[mRepo.readPos];
    chunk = mRepo.packBuffer[mRepo.readPos];
    //cerr << "read pos is " << mRepo.readPos << endl;
    //cerr << (char*)chunk->data.Pointer() << endl;

    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();

    double t = get_wall_time();
    //data format for from dsrc to fastp
    data->count = dsrc::fq::chunkFormat(chunk, data->data, true);

    config->costFormat += get_wall_time() - t;    //cerr << (char*)chunk->data.Pointer() << endl;
    fastqPool->Release(chunk);

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processSingleEnd(data, config);


}

void SingleEndProcessor::producerTask() {
    double t0 = get_wall_time();
    if (mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    //Read** data = new Read*[PACK_SIZE];
    //memset(data, 0, sizeof(Read*)*PACK_SIZE);
    //FastqReader reader(mOptions->in1, true, mOptions->phred64);
    //dsrc::fq::FastqDataPool* fastqPool = new dsrc::fq::FastqDataPool(32,1<<22);
    dsrc::fq::FastqFileReader *fileReader = new dsrc::fq::FastqFileReader(mOptions->in1);
    dsrc::fq::FastqReader *dataReader = new dsrc::fq::FastqReader(*fileReader, *fastqPool);


    dsrc::fq::FastqDataChunk *chunk;
    while ((chunk = dataReader->readNextChunk()) != NULL) {
        //cerr << "chunk content in producer =======================" << endl;
        //cerr << (char*) chunk->data.Pointer() << endl;
        producePack(chunk);
        //cerr << (char*)mRepo.packBuffer[0]->data.Pointer() << endl;
        //usleep(200);
        while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
            //cerr<<"sleep"<<endl;
            slept++;
            usleep(100);
        }

    }



    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if (mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();
    delete fileReader;
    delete dataReader;
    printf("producer cost %.5f\n", get_wall_time() - t0);
}

void SingleEndProcessor::consumerTask(ThreadConfig *config) {
//    cost1 = 0;
//    cost2 = 0;
//    cost3 = 0;
//    cost4 = 0;
//    cost5 = 0;
//    cost6 = 0;
//    cost7 = 0;
//    cost8 = 0;
//    cost9 = 0;
//    cost10 = 0;
//    cost11 = 0;
//    cost11b = 0;
//    cost12 = 0;
//    totCnt = 0;
    while (true) {
        if (config->canBeStopped()) {
            mFinishedThreads++;
            break;
        }
        while (mRepo.writePos <= mRepo.readPos) {
            if (mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
            mFinishedThreads++;
            if (mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if (mProduceFinished) {
            if (mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " +
                             to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        } //--[haoz:] I think it 3GS can add here
    }

//
//    printf("total getPreStats1()->statRead(or1) ====: %.5f\n", cost1);
//    printf("total mDuplicate->statRead(or1) ========: %.5f\n", cost2);
//    printf("total mOptions->indexFilter()  =========: %.5f\n", cost3);
//    printf("total mUmiProcessor->process(or1) ======: %.5f\n", cost4);
//    printf("total mFilter->trimAndCut() ============: %.5f\n", cost5);
//    printf("total PolyX::trimPolyG() ===============: %.5f\n", cost4);
//    printf("total trimBySequence ===================: %.5f\n", cost7);
//    printf("total r1->resize() =====================: %.5f\n", cost8);
//    printf("total mFilter->passFilter(r1) ==========: %.5f\n", cost9);
//    printf("total addFilterResult(result) ==========: %.5f\n", cost10);
//    printf("total outstr += r1->toString() =========: %.5f\n", cost11b);
//    printf("total getPostStats1()->statRead(r1) ====: %.5f\n", cost11);
//    printf("total delete r1 ========================: %.5f\n", cost12);
//    printf("total statReadCost =====================: %.5f\n", cost11 + cost11);
//    printf("total costTotel ========================: %.5f\n",
//           cost1 + cost2 + cost3 + cost4 + cost5 + cost6 + cost7 + cost8 + cost9 + cost10 + cost11b + cost11 + cost12);
//    printf("total cost =============================: %.5f\n", cost);
//    printf("total  =================================: %d\n", totCnt);

//TODO check another !
//    fstream out;
//    fstream in;
//
//    string outFileName1 = "oneStatsCheckOfThread" + to_string(config->getThreadId());
//    string outFileName2 = "kMerCheckOfThread" + to_string(config->getThreadId());
//    string outFileName3 = "mCounts" + to_string(config->getThreadId());
//    string outFileName4 = "mGC" + to_string(config->getThreadId());
//
//    long *tmpKmer = config->getPreStats1()->mKmer;
//    int tmpKmerSize = config->getPreStats1()->mKmerBufLen;
//    long *tmpStats = config->getPreStats1()->mCycleTotalQual;
//    int tmpStatsSize = config->getPreStats1()->mBufLen;
//    uint16 *tmpCounts = mDuplicate->mCounts;
//    int tmpCountSize = mDuplicate->mKeyLenInBit;
//    uint8 *tmpGC = mDuplicate->mGC;
//    int tmpGCSize = mDuplicate->mKeyLenInBit;
//
//    printf("Start checking ...");
//    printf("Open file : %s\n", outFileName1.c_str());
//
//    // Read sampling file
//    long *temp = new long[tmpStatsSize];
//    in.open(outFileName1.c_str(), ios::in | ios::binary);
//    if (!in) {
//        printf("Can't open file \"%s\"\n", outFileName1.c_str());
//    } else {
//        in.seekg(0, ios::beg);
//        in.read(reinterpret_cast<char *>(temp), tmpStatsSize * sizeof(long));
//        printf("=================================================\n");
//        for (int i = 0; i < tmpStatsSize; i++) {
//            if (temp[i] != tmpStats[i]) {
//                printf("GG on test %d  STD : %ld   Now : %ld\n", i, temp[i], tmpStats[i]);
//            }
////            printf("%ld\n", temp[i]);
//        }
//        printf("=================================================\n");
//    }
//
//
//    printf("=================================================\n");
//    for (int i = 0; i < tmpStatsSize; i++) {
//        printf("%ld ", tmpStats[i]);
//        if ((i + 1) % 20 == 0)printf("\n");
//    }
//    printf("=================================================\n");
//
//    out.open(outFileName1.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpStats), tmpStatsSize * sizeof(long));
//    out.close();
//
//    out.open(outFileName2.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpKmer), tmpKmerSize * sizeof(long));
//    out.close();
//
//    out.open(outFileName3.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpCounts), tmpCountSize * sizeof(long));
//    out.close();
//
//    out.open(outFileName4.c_str(), ios::out | ios::binary);
//    out.seekp(0, ios::beg);
//    out.write(reinterpret_cast<char *>(tmpGC), tmpGCSize * sizeof(long));
//    out.close();


    if (mFinishedThreads == mOptions->thread) {
        if (mLeftWriter)
            mLeftWriter->setInputCompleted();
    }

    if (mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writeTask(WriterThread *config) {
    double t0 = get_wall_time();
    while (true) {
        if (config->isCompleted()) {
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if (mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
    printf("write cost %.5f\n", get_wall_time() - t0);
}
