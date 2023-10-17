#pragma once // This is a common directive to prevent multiple inclusions of the header file

#include <vector>
#include <queue>

using namespace std;

class multiLevelQueue
{
private:
    vector<queue<int>> multi_queue;
    int max_level;

public:
    multiLevelQueue();
    ~multiLevelQueue();
    void initial(int n);
    void push(int i, int level);
    int get_next();
    bool empty();
};
