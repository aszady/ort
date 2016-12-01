#include <includes/header.hpp>

template<class V, class Trans>
struct GSegTree {
  // Assumes that:
  // Mix(a, Mix(b, c)) == Mix(Mix(a, b), c)
  // Does not assume:
  // Mix(a, b) == Mix(b, a)
  
  // V is default constructible
  // V is movable
  // V is copy-constructible
  GSegTree(vector<V> initial, function<V(const V&, const V&)> mix)
  : S(initial.size()), Mix(mix) {  assert(!initial.empty()); assert(mix);
    static_assert(is_default_constructible<V>::value, "V def-con");
    SR = 1 << (__builtin_clz(1) - __builtin_clz((S - 1) | 1) + 1);
    D = vector<V>(SR*2);  // Assuming V is default-constructible
    T = vector<Trans>(SR, Trans::neutral());    
    
    move(initial.begin(), initial.end(), D.begin() + SR);
    for (int i = SR; --i>0;) 
      if (valid(i))  // This enforces neutral-correctness
        D[i] = Mix(D[i*2], D[i*2+1]);
  } 
  
  GSegTree(const GSegTree&) = delete; 
  GSegTree& operator=(const GSegTree&) = delete;
  
  GSegTree(const GSegTree& other, function<V(const V&, const V&)> mix)
  : S(other.S)
  , SR(other.SR)
  , D(other.D)
  , T(other.T)
  , Mix(mix)
  {
    assert(mix);
  }
  
  int size() const { return S; }
  
  void assertValid() const { assert(size()>0); assert(Mix); }
 
  V query(int a, int b) const {  assert(a>=0); assert(a<S);
                                 assert(b>=0); assert(b<S);
                                 assert(a<=b);
    assert(Mix);
    const auto cq = bases(a, b).cq;
    assert(!cq.empty());
    auto it = cq.begin();
    assert(*it >= SR);
    V v(D[*(it++)]);
    
    for (; it != cq.end(); ++it){ 
      if (*it < SR) {
        v = Mix(v, T[*it].combine(D[*it], rl(*it)));
      } else {        
        v = Mix(v, D[*it]);
      }
    }
    return v;
  }
  
  void queryCustom(int a, int b, const function<void(const V&, int /*ra*/, int /*rb*/)>& fn) const {
    assert(Mix);
    for (int x : bases(a, b).cq)
      fn(D[x], ra(x), rb(x));
  }
  
  void apply(int a, int b, Trans trans) {  assert(a>=0); assert(a<S);
                                           assert(b>=0); assert(b<S);
                                           assert(a<=b);       
    auto ba = bases(a, b);    
    for (int x : ba.cq) {
      if (x < SR) 
        T[x] = trans.move(ra(x)-a, rl(x)-(b-a+1)).compose(T[x]);
      else 
        D[x] = trans.apply(D[x], ra(x)-a);
    }
    
    for (int x : ba.pq) {
      if (x < SR/2)
        D[x] = Mix(T[x*2  ].combine(D[x*2  ], rl(x*2  )),
                   T[x*2+1].combine(D[x*2+1], rl(x*2+1)));
      else  // if (x < SR)
        D[x] = Mix(D[x*2], D[x*2+1]);
    }
  }
  
  pair<typename vector<V>::const_iterator, typename vector<V>::const_iterator> getAll() const {
    return {D.begin() + SR, D.begin() + SR + S};
  }
  
  void resetMix(function<V(const V&, const V&)> mix) {
    Mix = mix;
  }
  
 private:
  struct B { vector<int> pq, cq; };
  
  B bases(int a, int b) const {  assert(a>=0); assert(a<=b); assert(b<S);
    B ba;
    auto& pq = ba.pq;
    auto& cq = ba.cq;  // In proper order (inorder)
    vector<int> rcq;    

    int u = a+SR, v = b+SR;
    cq.push_back(u);
    if (u != v) 
      rcq.push_back(v);
      
    while (u/2 != v/2) {
      if (!(u&1)) cq.push_back(u+1); 
      if (v&1) rcq.push_back(v-1); 
      pq.push_back(u/=2);
      if (valid(v/=2)) pq.push_back(v);
    }
    while (u/=2) {
      if (!valid(u)) break;
      pq.push_back(u);
    }

    cq.insert(cq.end(), rcq.rbegin(), rcq.rend());

    // Propagate Trans down and make them all TN
    bool changed = false;
    for (auto it = pq.rbegin(); it != pq.rend(); ++it) {
      int i = *it;
      if (T[i].isNeutral()) continue;
      if (i < SR/2) {
        T[i].move(0, -rl(i*2+1)).compose(&T[i*2]);
        T[i].move(rl(i*2), -rl(i*2)).compose(&T[i*2+1]);
      } else if (i < SR) {
        T[i].apply(&D[i*2  ], 0);
        T[i].apply(&D[i*2+1], 1);
      }
      changed = true;
      T[i] = Trans::neutral();  // Assuming i < SR
    }

    // Recalculate D upwards
    if (changed)
      for (int x : pq)
        if (x < SR)
          D[x] = Mix(D[x*2], D[x*2+1]);
    
    return ba;
  }

  int rl(int x) const {  assert(x>0 && x < 2*SR);  // Length of an interval pointed by node x
    return 1<<(__builtin_clz(x)-__builtin_clz(SR));
  }
  int ra(int x) const {  assert(x>0 && x < 2*SR);  // Start idx of an interval pointed by node x
    return x*rl(x)-SR;
  }
  int rb(int x) const {  assert(x>0 && x < 2*SR);  // End idx of an interval pointed by node x
    return ra(x)+rl(x)-1;
  }   
  bool valid(int x) const {  assert(x>0 && x < 2*SR);
    return rb(x) < S;
  }

  int S, SR;
  mutable vector<V> D;
  mutable vector<Trans> T;
  function<V(const V&, const V&)> Mix;
};
