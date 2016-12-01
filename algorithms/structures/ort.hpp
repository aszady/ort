#include "algorithms/structures/gsegtree.hpp"

template<class V>
using MixT = function<V(V, V)>;


struct Presorted {};  

template<size_t Dim, size_t IthDim, class V, class GSegTreeV, class GSegTreeTrans>
class ORTStructTraits {  
 public:
  int size() const { assert(segTree_); return segTree_->size(); }  
  
 protected:   
  
  ORTStructTraits(const std::vector<GSegTreeV>& initial,
                  const MixT<V>& mix,
                  function<GSegTreeV(const GSegTreeV&, const GSegTreeV&)> gMix)
  : Mix(mix)
  , segTree_(new GSegTree<GSegTreeV, GSegTreeTrans>(initial, gMix))
  {}
  
  ORTStructTraits(const ORTStructTraits& other,  
    const MixT<V>& mix,
    function<GSegTreeV(const GSegTreeV&, const GSegTreeV&)> gMix)
  : Mix(mix)
  , segTree_(other.segTree_ ? new GSegTree<GSegTreeV, GSegTreeTrans>(*other.segTree_, gMix) : nullptr)
  {}
  
  ORTStructTraits() = default;
  
  void swap(ORTStructTraits&& other,
            function<GSegTreeV(const GSegTreeV&, const GSegTreeV&)> gMix,
            function<GSegTreeV(const GSegTreeV&, const GSegTreeV&)> gMixOther) {    
    std::swap(Mix, other.Mix);
    std::swap(segTree_, other.segTree_);
    if (segTree_) {
      assert(gMix);
      segTree_->resetMix(gMix);
    }
    if (other.segTree_) {
      assert(gMixOther);
      other.segTree_->resetMix(gMixOther);
    }
  }
  
  template<class Locator>
  pair<int, int> locate(const V& a, const V& b) const {
    auto range = segTree_->getAll();
    auto ita = lower_bound(range.first, range.second, a, Locator{});
    auto itb = lower_bound(range.first, range.second, b, Locator{});
    int pa = ita - range.first;
    int pb = itb - range.first;
    
    return {pa, pb};
  }  
 
  MixT<V> Mix;
  unique_ptr<GSegTree<GSegTreeV, GSegTreeTrans>> segTree_;
};

template<size_t Dim, size_t IthDim, class V, class DimCmp, class Trans>
class ORTStruct;


template<class E>
struct EmptyTrans {  
  void apply(E* /*ort*/, int /*dx*/) { assert(false); }
  E combine(const E& e, int /*len*/) { assert(false); return e; }  
  void compose(EmptyTrans* /*t*/) { assert(false); }
  EmptyTrans move(int /*dx*/, int /*dl*/) { assert(false); return {}; }
  static EmptyTrans neutral() { return {}; }
  bool isNeutral() const { return true; }
};

template<size_t Dim, size_t IthDim, class V, class DimCmp, class Trans>
class ORTStruct : public ORTStructTraits<
    Dim, IthDim, V,
    ORTStruct<Dim, IthDim-1, V, DimCmp, Trans>,
    EmptyTrans<ORTStruct<Dim, IthDim-1, V, DimCmp, Trans>>
> { 
  using Base = ORTStructTraits<
    Dim, IthDim, V,
    ORTStruct<Dim, IthDim-1, V, DimCmp, Trans>,
     EmptyTrans<ORTStruct<Dim, IthDim-1, V, DimCmp, Trans>>
  >;
  using NextORT = ORTStruct<Dim, IthDim-1, V, DimCmp, Trans> ;
  
 public:
    
  explicit ORTStruct(const vector<V>& initial, Presorted, const MixT<V>& mix)
  : Base(
    intoSingleOrtStructs(initial, mix),
    mix,
    bind(&ORTStruct::mixer, this, std::placeholders::_1, std::placeholders::_2)) {  assert(!initial.empty()); assert(mix);
    }
    
  
  explicit ORTStruct(const vector<V>& initial, const MixT<V>& mix)
  : ORTStruct(sorted(initial), Presorted{}, mix)
  {}
  
  ORTStruct() = default; // To be default-constructible by GSegTree. GSegTree guarantees that this object won't be used
  
  ORTStruct(const ORTStruct& other)
  : Base(
      other,
      other.Mix,
      bind(&ORTStruct::mixer, this, std::placeholders::_1, std::placeholders::_2))
  {}
  
  ORTStruct& operator=(ORTStruct&& other) {
    Base::swap(std::move(other),
      bind(&ORTStruct::mixer, this, std::placeholders::_1, std::placeholders::_2),
      bind(&ORTStruct::mixer, &other, std::placeholders::_1, std::placeholders::_2));
    return *this;
  }
  
  void assertValid() const {
     assert(Base::size()>0);
     assert(Base::Mix);
     assert(Base::segTree_);
     assert(Base::size() == Base::segTree_->size());
     Base::segTree_->assertValid();
  }
  
  struct QueryLocatorComparator {
    bool operator()(const NextORT& o, const V& v) const {
      o.assertValid();
      assert(o.size() == 1);
      return DimCmp{}.template precedes<IthDim>(o.querySingleton(), v);
    }
  };
  
  template<class Debugger>
  V query(const V& a, const V& b, bool& any, Debugger& debugger) const {  assert(Base::segTree_);    
    any = false;      
    auto range = Base::template locate<QueryLocatorComparator>(a, b);
    int pa = range.first;
    int pb = range.second - 1;  
  
    debugger.onQueryStart(IthDim, a, b);
    if (pa > pb) return V{};
    
    V v;
    Base::segTree_->queryCustom(pa, pb, [&, this](const NextORT& o, int ra, int rb) {
      debugger.onPerspectiveSet(IthDim, 
        Base::segTree_->query(ra, ra).querySingleton(), Base::segTree_->query(rb, rb).querySingleton()
      );
      bool any_rec = false;
      V rv = o.query(a, b, any_rec, debugger);
      if (any_rec) {
        v = !any ? any=true, rv : Base::Mix(std::move(v), std::move(rv));
      }
    });
    return v;
  }
  
  void apply(const V& a, const V& b, const Trans& t) { assert(Base::segTree_);
    Base::template locateAndQuery<QueryLocatorComparator>(a, b, 
      [&, this](int pa, int pb) {
        Base::segTree_->queryCustom(pa, pb, [&, this](const NextORT& o, int ra, int rb) {
          o.apply(a, b, t);
        });
      }
    );
  }
  
  V querySingleton() const { assert(Base::segTree_ && Base::Mix); assert(Base::size() == 1);
    return Base::segTree_->query(0, 0).querySingleton();
  }
  
  vector<V> getAll() const { assert(Base::segTree_ && Base::Mix);
    auto range = Base::segTree_->getAll();
    vector<V> sum;
    for (auto it = range.first; it != range.second; ++it) {
      const NextORT& no = *it;   
      assert(no.size() == 1);
      sum.push_back(no.querySingleton());
    }
    return sum;    
  }
  
 private:
  static vector<V> sorted(vector<V> v) {    
    sort(v.begin(), v.end(), [](const V& v1, const V& v2) { return DimCmp{}.template precedes<IthDim>(v1, v2); });
    return v;
  }
  
  static vector<NextORT> intoSingleOrtStructs(const vector<V>& sorted, const MixT<V>& mix) {        
    vector<NextORT> result;
    for (const V& v : sorted) {
      result.push_back(NextORT({v}, mix));
      result.back().assertValid();
    }
    return result;
  }
  
  NextORT mixer(
    const NextORT& left,
    const NextORT& right) const 
  {
    vector<V> sum;
    auto l = left.getAll();
    auto r = right.getAll();
    merge(l.begin(), l.end(), r.begin(), r.end(), back_inserter(sum),
      [](const V& v1, const V& v2) { return DimCmp{}.template precedes<IthDim-1>(v1, v2); }
    );
    return NextORT{sum, Presorted{}, Base::Mix};
  }
};

template<size_t Dim, class V, class DimCmp, class Trans>
class ORTStruct<Dim, 0, V, DimCmp, Trans> : public ORTStructTraits<
    Dim, 0, V, V, Trans> {  
  using Base = ORTStructTraits<Dim, 0, V, V, Trans>;
 public:
  explicit ORTStruct(const vector<V>& initial, Presorted, const MixT<V>& mix)
  : Base(initial, mix, mix) {  assert(!initial.empty());  assert(mix);
  }
  
  explicit ORTStruct(const vector<V>& initial, const MixT<V>& mix)
  : ORTStruct(sorted(initial), Presorted{}, mix) {  assert(!initial.empty());  assert(mix);
  }
  
  ORTStruct() = default;  // To be default-constructible by GSegTree. GSegTree guarantees that this object won't be used
  
  ORTStruct(const ORTStruct& other)
  : Base(other, other.Mix, other.Mix)
  {}
  
  ORTStruct& operator=(ORTStruct&& other) {
    Base::swap(std::move(other), other.Mix, Base::Mix);
    return *this;
  }
  
  void assertValid() const {
    assert(Base::size()>0);
    assert(Base::Mix);
    assert(Base::segTree_);
    assert(Base::size() == Base::segTree_->size());
    Base::segTree_->assertValid();
  }
  
  struct QueryLocatorComparator {
    bool operator()(const V& v1, const V& v2) const {
      return DimCmp{}.template precedes<0>(v1, v2);
    }
  };
  
  template<class Debugger>
  V query(const V& a, const V& b, bool& any, Debugger& debugger) const { assert(Base::segTree_);    
    any = false;
    const auto range = Base::template locate<QueryLocatorComparator>(a, b);
    int pa = range.first;
    int pb = range.second - 1;
    
    debugger.onQueryStart(0, a, b);
    
    if (pa>pb) return V{};
    
    V v;
    Base::segTree_->queryCustom(pa, pb, [&, this](const V& val, int ra, int rb) {
      debugger.onPerspectiveSet(0, 
        Base::segTree_->query(ra, ra), Base::segTree_->query(rb, rb)
      );
      debugger.onLastDimFound(val);
      v = !any ? any=true, val : Base::Mix(std::move(v), val);
    });
    return v;
  }  
  
  V querySingleton() const { assert(Base::segTree_ && Base::Mix); assert(Base::size() == 1);
    return Base::segTree_->query(0, 0);
  }
  
  vector<V> getAll() const { assert(Base::segTree_ && Base::Mix);
    auto range = Base::segTree_->getAll();
    return {range.first, range.second};
  }
  
  void apply(const V& a, const V& b, const Trans& t) const { assert(Base::segTree_); 
    return Base::template locateAndQuery<QueryLocatorComparator>(a, b,
      [this](int pa, int pb) {
        return Base::segTree_->query(pa, pb);
      }
    );
  }  
  
  void applyAll(const Trans& t) { assert(Base::segTree_ && Base::Mix);
    return Base::segTree_->apply(0, Base::segTree_->size()-1, t);
  }
  
 private:
  static vector<V> sorted(vector<V> s) {
    sort(s.begin(), s.end(), [](const V& v1, const V& v2) { return DimCmp{}.template precedes<0>(v1, v2); });
    return s;
  }
};

template<class V>
struct ORTEmptyDebugger {
  void onQueryStart(size_t /*dim*/, const V& /*v1*/, const V& /*v2*/) {}
  void onPerspectiveSet(size_t /*dim*/, const V& /*v1*/, const V& /*v2*/) {}
  void onLastDimFound(const V&) {}
};

template<size_t Dim, class V, class DimCmp, class Trans>
class ORT {
 public:  
  explicit ORT(const vector<V>& initial, const function<V(V, V)>& mix)
  : struct_(initial, mix) {}
  
  template<class Debugger>
  V query(const V& a, const V& b, bool& any, Debugger& debugger) const {
    return struct_.query(a, b, any, debugger);
  }
  
  V query(const V& a, const V& b, bool& any) const {
    ORTEmptyDebugger<V> debugger;
    return query(a, b, any, debugger);
  }
  
  vector<V> getAll() const {
    return struct_.getAll();
  }
  
 private:
  ORTStruct<Dim, Dim-1, V, DimCmp, Trans> struct_;
};
