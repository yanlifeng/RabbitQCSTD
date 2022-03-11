#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>

using namespace std;

struct node {
    string s1, s2;

    bool operator<(const node A) const {
        return s1 < A.s1;
    }
};

int main(int argc, char *argv[]) {

    auto in_name = argv[1];
    auto out_name = argv[2];


    string s1, s2, s3 = "";
    vector<node> A;

    freopen(in_name, "r", stdin);
    while (cin >> s1 >> s2) {
        node t;
        t.s1 = s1;
        t.s2 = s2;
        A.push_back(t);
    }

    fclose(stdin);
    sort(A.begin(), A.end());
    freopen(out_name, "w", stdout);
    for (auto item:A) {
        cout << item.s1 << " " << item.s2 << endl;
    }
    fclose(stdout);
    return 0;
}