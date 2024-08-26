#ifndef __DSU_H__
#define __DSU_H__

#include <unordered_map>
#include <cassert>

template <typename T>
class DSU
{
public:
    void make_set(T t)
    {
        if (_parent.find(t) != _parent.end())
            return;
        _parent[t] = t;
        _rank[t] = 0;
    }

    T find_set(T t)
    {
        if (_parent[t] == t)
            return t;
        return _parent[t] = find_set(_parent[t]);
    }

    void union_set(T t1, T t2)
    {
        t1 = find_set(t1);
        t2 = find_set(t2);
        if (t1 != t2)
        {
            if (_rank[t1] < _rank[t2])
                std::swap(t1, t2);
            _parent[t2] = t1;
            if (_rank[t1] == _rank[t2])
                ++_rank[t1];
        }
    }

private:
    std::unordered_map<T, T> _parent;
    std::unordered_map<T, unsigned int> _rank;
};

#endif