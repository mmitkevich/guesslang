#ifndef QQFUNCTIONAL_H
#define QQFUNCTIONAL_H
/*
template<typename K, typename T, typename F, typename V>
inline V reduce(QMap<K,T>& dict, const V& start, F reduce) {
    for(auto it = dict.cbegin(); it!=dict.cend(); it++) {
        start = reduce(it.key(), it.value(), start);
    }
}

template<typename K, typename T, typename F>
inline typename QHash<K,T>::iterator loop(QHash<K,T>& dict, F f) {
    for(auto it = dict.begin(); it!=dict.end(); it++) {
        if(!f(it.key(), it.value()))
            return it;
    }
}

template<typename K, typename T, typename F>
inline typename QHash<K,T>::const_iterator loop(const QHash<K,T>& dict, F f) {
    for(auto it = dict.cbegin(); it!=dict.cend(); it++) {
        if(f(it.key(), it.value()))
            return it;
    }
}


template<typename K, typename T, typename F>
inline typename QMap<K,T>::iterator loop(QMap<K,T>& dict, F f) {
    for(auto it = dict.begin(); it!=dict.end(); it++) {
        if(f(it.key(), it.value()))
            return it;
    }
}

template<typename K, typename T, typename F>
inline typename QMap<K,T>::const_iterator loop(const QMap<K,T>& dict, F f) {
    for(auto it = dict.cbegin(); it!=dict.cend(); it++) {
        if(!f(it.key(), it.value()))
            return it;
    }
}

template<typename V, typename T, typename F>
QPair<int, V> max_element(const QList<T>& list, F f, const V& minval = std::numeric_limits<V>::min, bool first=true) {
    QPair<int,V> best;
    best.first = -1;
    best.second = minval;
    for(int i = first ? 0 : list.size()-1; first ? i<list.size() : i>=0; first?i++:i--) {
        auto val = f(i, list[i]);
        if(best.second < val) {
            best.second = val;
            best.first = i;
        }
    }
    return best;
}

// transform K->T into V->R and findout R with maximum V
template<typename V, typename R, typename K, typename T, typename F>
QPair<V, R> max(const QMap<K, T>& dict, F f, QPair<V,R> &&min = QPair<V,R>(), bool first=true) {
    QPair<V,R> best = std::move(min);
    for(auto it = dict.cbegin(); it!=dict.cend(); it++) {
        auto pair = f(it.key(), it.value()); // qpair: key compares
        if(first && best.first < pair.first || !first && best.first <= pair.first) {
            best = std::move(pair);
        }
    }
    return best;
}

// reduce {K->T} into R
template<typename R, typename K, typename T, typename F>
R reduce(const QMap<K, T>& dict, F f, R &&start = R(), bool first=true) {
    R result = start;
    for(auto it = dict.cbegin(); it!=dict.cend(); it++) {
        if(!f(it.key(), it.value(), result))
            return result;
    }
    return result;
}

#define ldef(expr, ...)  [&](__VA_ARGS__){ return expr; }
*/
#endif // QQFUNCTIONAL_H
