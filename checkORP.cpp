#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>

using namespace std;

struct node {
    string s1, s2;

    bool operator<(const node A) const {
        if (s1 == A.s1)return s2 < A.s2;
        return s1 < A.s1;
    }
};

int main() {

    string s1, s2;
    vector<node> A;
    vector<node> B;
    int ok;

    freopen("ORP.log", "r", stdin);
    while (cin >> s1) {
        node t;
        t.s1 = s1;
        A.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    freopen("../RabbitQCPlus/ORP.log", "r", stdin);
    while (cin >> s1) {
        node t;
        t.s1 = s1;
        B.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    if (A.size() != B.size()) {
        printf("GG on size %d %d\n", A.size(), B.size());
//        exit(0);
    }
    ok = 1;
    for (int i = 0; i < min(A.size(), B.size()); i++) {
        if (A[i].s1 != B[i].s1)ok = 0;
    }
    if (ok == 0)printf("GG in check ORP\n");
    else printf("AC in check ORP\n");

    A.clear();
    B.clear();

    freopen("ORP2.log", "r", stdin);
    while (cin >> s1 >> s2) {
        node t;
        t.s1 = s1;
        t.s2 = s2;
        A.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    freopen("../RabbitQCPlus/ORP2.log", "r", stdin);
    while (cin >> s1 >> s2) {
        node t;
        t.s1 = s1;
        t.s2 = s2;
        B.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    if (A.size() != B.size()) {
        printf("GG on size %d %d\n", A.size(), B.size());
//        exit(0);
    }
    ok = 1;
    for (int i = 0; i < min(A.size(), B.size()); i++) {
        if (A[i].s1 != B[i].s1) {
            ok = 0;
//            cout << A[i].s1 << " " << B[i].s1 << endl;
        }
        if (A[i].s2 != B[i].s2) {
            ok = 0;
//            cout << A[i].s2 << " " << B[i].s2 << endl;

        }
    }
    if (ok == 0)printf("GG in check ORP2\n");
    else printf("AC in check ORP2\n");
    fclose(stdin);
    cin.clear();


    A.clear();
    B.clear();

    freopen("ORP3.log", "r", stdin);
    while (cin >> s1 >> s2) {
        node t;
        t.s1 = s1;
        t.s2 = s2;
        A.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    freopen("../RabbitQCPlus/ORP3.log", "r", stdin);
    while (cin >> s1 >> s2) {
        node t;
        t.s1 = s1;
        t.s2 = s2;
        B.push_back(t);
    }
    fclose(stdin);
    cin.clear();
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    if (A.size() != B.size()) {
        printf("GG on size %d %d\n", A.size(), B.size());
//        exit(0);
    }
    ok = 1;
    for (int i = 0; i < min(A.size(), B.size()); i++) {
        if (A[i].s1 != B[i].s1) {
            ok = 0;
//            cout << A[i].s1 << " " << B[i].s1 << endl;
        }
        if (A[i].s2 != B[i].s2) {
            ok = 0;
//            cout << A[i].s2 << " " << B[i].s2 << endl;

        }
    }
    if (ok == 0)printf("GG in check ORP3\n");
    else printf("AC in check ORP3\n");
    fclose(stdin);

    return 0;
}
