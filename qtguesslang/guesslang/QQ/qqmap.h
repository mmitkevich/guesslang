#ifndef QQMAP_H
#define QQMAP_H

#include <QtCore/QMap>
#include <QtCore/QHash>

/// QMap<int,QString> a(1=>'a', 3=>'b'), b(1=>'c', 4=>'d');
/// a<<b
/// a == (1=>'c',3=>'b',4=>'d')
template<typename K, typename T>
inline QMap<K,T>& operator<<(QMap<K,T> &lval, const QMap<K,T>& rval) {
    for(auto i = rval.cbegin(); i!=rval.cend(); ++i) {
        lval[i.key()] =  i.value();
    }
    return lval;
}


template<typename K, typename T>
inline QHash<K,T>& operator<<(QHash<K,T> &lval, const QHash<K,T>& rval) {
    for(const auto & i = rval.cbegin(); i!=rval.cend(); i++) {
        lval.insert(i.key(), i.value());
    }
    return lval;
}


template<typename K, typename T>
inline QMap<K,T>& operator<<(QMap<K,T> &lval, const QPair<K,T> &&pair) {
    lval[pair.first] = std::move(pair.second);
    return lval;
}

template<typename K, typename T>
inline QHash<K,T>& operator<<(QHash<K,T> &lval, const QPair<K,T> &pair) {
    lval.insert(pair.first, pair.second)    ;
    return lval;
}


//template<typename K, typename T>
//inline QMap<K,T>& operator>>(QMap<K,T> &lval, T& rval) {
//    rval = lval.take(lval.firstKey());
//    return lval;
//}

template<class Map>
struct QQItems {
    typedef typename Map::const_iterator Iterator;
    typedef typename Map::mapped_type value_type;
    typedef typename Map::key_type key_type;

    const Map &map_;

    QQItems(const Map & map_) : map_(map_) {}

    struct iterator {
        Iterator mapIterator;
        iterator(const Iterator &mapIterator_): mapIterator(mapIterator_) {}
        Iterator operator*() {
            return mapIterator;
        }
        iterator & operator++() {
            ++mapIterator;
            return *this;
        }
        bool operator!=(const iterator & other) {
            return this->mapIterator != other.mapIterator;
        }
    };

    int size() const {
        return map_.size();
    }

    const value_type &operator[](const key_type& key) const {
        return map_[key];
    }

    iterator begin() const {
        return map_.cbegin();
    }

    iterator end() const {
        return map_.cend();
    }

    iterator cbegin() const {
        return map_.cbegin();
    }

    iterator cend() const {
        return map_.cend();
    }
};

/// for(auto it = items(myQmap)) cout << it.key() << "->" << it.value()
template<class Map> QQItems<Map> qqitems(Map & map) {
    return QQItems<Map>(map);
}


#endif // QQMAP_H
