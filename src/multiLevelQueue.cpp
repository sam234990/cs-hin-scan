#include "multiLevelQueue.h"

multiLevelQueue::multiLevelQueue()
{
}

multiLevelQueue::~multiLevelQueue()
{
}

void multiLevelQueue::initial(int n)
{
    multi_queue.resize(n);
    max_level = 0;
    for (auto &levelQueue : multi_queue)
        levelQueue = std::queue<int>();
}

void multiLevelQueue::push(int i, int level)
{
    if (level > max_level)
    {
        max_level = level;
    }
    multi_queue[level].push(i);
}

int multiLevelQueue::get_next()
{
    while (max_level >= 0)
    {
        if (multi_queue[max_level].empty())
        {
            max_level--;
            continue;
        }
        else
        {
            int res_vertex_id = multi_queue[max_level].front();
            multi_queue[max_level].pop();
            return res_vertex_id;
        }
    }
    if (max_level <= 0 && multi_queue[0].empty())
    {
        return -1;
    }
    return -1;
}

bool multiLevelQueue::empty()
{
    if (max_level > 0)
        return false;
    else
    {
        if (multi_queue[0].empty())
            return true;
        else
            return false;
    }
}
