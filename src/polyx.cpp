#include "polyx.h"

PolyX::PolyX() {
}


PolyX::~PolyX() {
}

void PolyX::trimPolyG(Read *r1, Read *r2, FilterResult *fr, int compareReq) {
    trimPolyG(r1, fr, compareReq);
    trimPolyG(r2, fr, compareReq);
}

void PolyX::trimPolyG(Read *r, FilterResult *fr, int compareReq) {
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = r->mSeq.mStr.c_str();

    int rlen = r->length();

    int mismatch = 0;
    int i = 0;
    int firstGPos = rlen - 1;
    for (i = 0; i < rlen; i++) {
        if (data[rlen - i - 1] != 'G') {
            mismatch++;
        } else {
            firstGPos = rlen - i - 1;
        }

        int allowedMismatch = (i + 1) / allowOneMismatchForEach;
        if (mismatch > maxMismatch || (mismatch > allowedMismatch && i >= compareReq - 1))
            break;
    }

    if (i >= compareReq) {
        r->resize(firstGPos);
    }
}

void PolyX::trimPolyX(Read *r1, Read *r2, FilterResult *fr, int compareReq) {
    trimPolyX(r1, fr, compareReq);
    trimPolyX(r2, fr, compareReq);
}

void PolyX::trimPolyX(Read *r, FilterResult *fr, int compareReq) {
    static bool geok = 0;
    if (geok == 0) {
        geok = 1;
        printf("compareReq %d\n", compareReq);
    }
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char *data = r->mSeq.mStr.c_str();

    int rlen = r->length();

    int atcgNumbers[4] = {0, 0, 0, 0};
    char atcgBases[4] = {'A', 'T', 'C', 'G'};
    int pos = 0;
    for (pos = 0; pos < rlen; pos++) {
        switch (data[rlen - pos - 1]) {
            case 'A':
                atcgNumbers[0]++;
                break;
            case 'T':
                atcgNumbers[1]++;
                break;
            case 'C':
                atcgNumbers[2]++;
                break;
            case 'G':
                atcgNumbers[3]++;
                break;
            case 'N':
                atcgNumbers[0]++;
                atcgNumbers[1]++;
                atcgNumbers[2]++;
                atcgNumbers[3]++;
                break;
            default:
                break;
        }

        int cmp = (pos + 1);
        int allowedMismatch = min(maxMismatch, cmp / allowOneMismatchForEach);

        bool needToBreak = true;
        for (int b = 0; b < 4; b++) {
            if (cmp - atcgNumbers[b] <= allowedMismatch)
                needToBreak = false;
        }
        if (needToBreak && (pos >= allowOneMismatchForEach || pos + 1 >= compareReq - 1)) {
            break;
        }
    }

    // has polyX
    if (pos + 1 >= compareReq && pos < rlen) {
//        static int cnt = 0;
//        static bool okk = 0;
//        static std::ofstream out_stream_;
//        if (okk == 0) {
//            out_stream_.open("info.txt");
//            okk = 1;
//        }
//        out_stream_ << "name " << r->mName << " " << cnt << " " << pos << " ";
//        out_stream_ << r->mSeq.mStr << " ";
//        out_stream_ << atcgNumbers[0] << " " << atcgNumbers[1] << " " << atcgNumbers[2] << " " << atcgNumbers[3]
//                    << " ";
//        out_stream_ << rlen << " ||| ";
//
//        cnt++;
        // find the poly
        int poly;
        int maxCount = -1;
        for (int b = 0; b < 4; b++) {
            if (atcgNumbers[b] > maxCount) {
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
//        out_stream_ << poly << " " << maxCount << " ||| ";
        char polyBase = atcgBases[poly];
//        out_stream_ << polyBase << " ||| " << rlen - pos - 1 << " ||| ";
        while (data[rlen - pos - 1] != polyBase && pos >= 0)
            pos--;
//        out_stream_ << pos << std::endl;

        r->resize(rlen - pos - 1);
    }
}

bool PolyX::test() {
    Read r("@name",
           "ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT",
           "+",
           "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    PolyX::trimPolyX(&r, NULL, 10);
    r.print();
    return r.mSeq.mStr == "ATTTT";
}