#include "gogui.hpp"
#include <cstdio>
#include <chrono>
#include <fstream>
#include <iterator>
#include <queue>
#include <set>
#include <unordered_set>

#include "algorithms/structures/ort.hpp"

using gogui::Point;
using gogui::Line;
using namespace gogui::utils;

template<size_t Dim>
using NDPoint = std::array<double, Dim>;


template<size_t Dim>
struct DimCmpSingle {
  template<size_t IthDim>
  static bool precedes(const NDPoint<Dim>& p1, const NDPoint<Dim>& p2) {
    return p1[IthDim] < p2[IthDim];
  }
};

template<size_t Dim>
struct DimCmpInVec {
  template<size_t IthDim>
  static bool precedes(const vector<NDPoint<Dim>>& p1, const vector<NDPoint<Dim>>& p2) {
    assert(p1.size() == 1);
    assert(p2.size() == 1);
    return DimCmpSingle<Dim>{}.template precedes<IthDim>(p1.front(), p2.front());
  }
};

template<size_t Dim>
struct SumOfPointsMix {
  using V = vector<NDPoint<Dim>>;
  V operator()(V v1, V v2) {
    v1.reserve(v1.size() + v2.size());
    move(v2.begin(), v2.end(), back_inserter(v1));
    return v1;
  }
};

template<size_t Dim>
struct SumOfPointsTrans {
  vector<NDPoint<Dim>> combine(vector<NDPoint<Dim>> v, int /*len*/) { return v; }
  void apply(vector<NDPoint<Dim>>* /*v*/, int /*dx*/) { }
  void compose(SumOfPointsTrans* /*t*/) { }
  SumOfPointsTrans move(int /*dx*/, int /*dl*/) { return *this; }
  static SumOfPointsTrans neutral() { return SumOfPointsTrans{}; }
  bool isNeutral() const { return true; }
};

template<size_t Dim>
struct MaxValueMix {
  using V = NDPoint<Dim+1>;
  V operator()(V v1, V v2) {
    return v1[Dim] > v2[Dim] ? v1 : v2;
  }
};

template<size_t Dim>
struct MaxValueTrans {
  NDPoint<Dim+1>  combine(NDPoint<Dim+1> v, int /*len*/) { return v; }
  void apply(NDPoint<Dim+1> * /*v*/, int /*dx*/) { }
  void compose(MaxValueTrans* /*t*/) { }
  MaxValueTrans move(int /*dx*/, int /*dl*/) { return *this; }
  static MaxValueTrans neutral() { return MaxValueTrans{}; }
  bool isNeutral() const { return true; }
};

template<class V>
struct Debugger {
  void onQueryStart(size_t dim, int pa, int pb, const V& va, const V& vb) const {    
    cerr << "Querying at " << dim << ".d: [" << pa << " - " << pb << "] ";
    
    for (const auto& p : va) {
      cerr << "<";
      for (const auto& x : p) {
        cerr << x << ", ";
      }
      cerr << ">";
    }
    
    cerr << " - ";
    
    for (const auto& p : vb) {
      cerr << "<";
      for (const auto& x : p) {
        cerr << x << ", ";
      }
      cerr << ">";
    }
    cerr << endl;    
  }
  
  void onPerspectiveSet(size_t dim, int pa, int pb, const V& va, const V& vb) const {    
    cerr << "PS at " << dim << ".d: [" << pa << " - " << pb << "] ";
    
    for (const auto& p : va) {
      cerr << "<";
      for (const auto& x : p) {
        cerr << x << ", ";
      }
      cerr << ">";
    }
    
    cerr << " - ";
    
    for (const auto& p : vb) {
      cerr << "<";
      for (const auto& x : p) {
        cerr << x << ", ";
      }
      cerr << ">";
    }
    cerr << endl;    
  }
};

template<class T>
T timer(const std::string& label, const std::function<T()>& fn) {
  const auto start = std::chrono::system_clock::now();
  T v = fn();
  const auto end = std::chrono::system_clock::now();
  const std::chrono::duration<double> dt = end-start;
  std::cout << label << ": " << dt.count() << endl;
  return v;
}

template<>
void timer(const std::string& label, const std::function<void()>& fn) {
  timer(label, function<bool()>([fn] { fn(); return true; }));
}

template<size_t Dim, class DimCmp, class Mix, class Trans>
void testConstruction(
    const size_t n,
    const function<typename Mix::V()>& pointCreator,
    const function<typename Mix::V(array<double, Dim>)> queryPointCreator) {
  using V = typename Mix::V;
  using O = ORT<Dim, V, DimCmp, Trans>;   
  
  const auto data = timer("  Creating points", function<vector<V>()>([n, &pointCreator] {
    vector<V> data;
    for (size_t i = 0; i < n; ++i) {
      data.push_back(pointCreator());
    }
    return data;
  }));
  
  O tree = timer("  Constructing tree", function<O()>([&data] {
    return O(data, Mix{});
  }));
  
  auto queryFnFactory = [n, &tree, &queryPointCreator](double expected, size_t queries) -> function<void()> {
    const double dx = std::pow(expected/n, 1./Dim);
    return [dx, queries, &tree, &queryPointCreator]() mutable {
      while (queries--) {
        std::array<double, Dim> a, b;      
        for (size_t x = 0; x<Dim; ++x) {
          a[x] = (static_cast<double>(rand())/RAND_MAX) * (1-dx);
          b[x] = a[x] + dx;
        }      
        bool any;
        const V result = tree.query(queryPointCreator(a), queryPointCreator(b), any);
      }
    };
  };
  
  timer("  1000 small queries (expected 10 points)", queryFnFactory(10, 1));
  timer("  1000 med queries (expected 100 points)", queryFnFactory(100, 1));
  timer("  1000 big queries (expected 1000 points)", queryFnFactory(1000, 1));
}

template<size_t Dim, class DimCmp, class Mix, class Trans>
void testConstruction(
    const std::initializer_list<size_t> ns,
    const function<typename Mix::V()>& pointCreator,
    const function<typename Mix::V(array<double, Dim>)> queryPointCreator) {
  for (size_t n : ns) {
    std::cout << " N = " << n << std::endl;
    testConstruction<Dim, DimCmp, Mix, Trans>(n, pointCreator, queryPointCreator);
  }
}

template<size_t Dim>
struct RandomPointCreator {
  NDPoint<Dim> operator()() const {
    NDPoint<Dim> p;
    for (size_t x = 0; x<Dim; ++x)
      p[x] = static_cast<double>(rand())/RAND_MAX;
    return p;
  }
};

template<size_t Dim>
struct RandomPointCreatorInVec {
  vector<NDPoint<Dim>> operator()() const {
    return {RandomPointCreator<Dim>{}()};
  }
};

template<size_t Dim>
struct QueryPointCreatorValued {
  NDPoint<Dim+1> operator()(const std::array<double, Dim>& x) const {
    NDPoint<Dim+1> p;
    copy(x.begin(), x.end(), p.begin());
    return p;
  }
};

template<size_t Dim>
struct QueryPointCreatorInVec {
  vector<NDPoint<Dim>> operator()(const std::array<double, Dim>& x) const {
    return {x};
  }
};

template<size_t Dim>
void testConstructionForSumOfPoints(const std::initializer_list<size_t>& ns) {
  std::cout << "ORT " << Dim << "D with Mix = sum: " << std::endl;
  testConstruction<Dim, DimCmpInVec<Dim>, SumOfPointsMix<Dim>, SumOfPointsTrans<Dim>>(
    ns,
    RandomPointCreatorInVec<Dim>{},
    QueryPointCreatorInVec<Dim>{}
  ); 
}

template<size_t Dim>
void testConstructionForMax(const std::initializer_list<size_t>& ns) {
  std::cout << "ORT " << Dim << "D with Mix = max: " << std::endl;
  testConstruction<Dim, DimCmpSingle<Dim+1>, MaxValueMix<Dim>, MaxValueTrans<Dim>>(
    ns,
    RandomPointCreator<Dim+1>{},
    QueryPointCreatorValued<Dim>{}
  );
}

template<size_t Dim>
struct GoGuiVisualizer;

template<>
struct GoGuiVisualizer<2> {
  struct DrawnBox {
    DrawnBox(double xx1, double yy1, double xx2, double yy2, const string& c)
    : x1(xx1), y1(yy1), x2(xx2), y2(yy2)
    , l1({x1, y1}, {x2, y1})
    , l2({x1, y1}, {x1, y2})
    , l3({x1, y2}, {x2, y2})
    , l4({x2, y1}, {x2, y2})
    , al1{(l1.setColor(c), l1)}
    , al2{(l2.setColor(c), l2)}
    , al3{(l3.setColor(c), l3)}
    , al4{(l4.setColor(c), l4)}
    {}
    double x1, y1, x2, y2;
    gogui::Line l1, l2, l3, l4;
    gogui::ActiveLine al1, al2, al3, al4;
  };
  
  void onQueryStart(size_t dim, const vector<NDPoint<2>>& va, const vector<NDPoint<2>>& vb) {
    const auto& a = va.front();
    const auto& b = vb.front();
    
    if (dim == 1) {
      DrawnBox box(0, a[1], 1, b[1], "red");      
      gogui::snapshot("Dim=1/2");
    } else if (dim == 0) {
      DrawnBox box(a[0], boxes_[1]->y1, b[0], boxes_[1]->y2, "red");
      gogui::snapshot("Dim=0/2");
    }
  }
  
  void onPerspectiveSet(size_t dim, const vector<NDPoint<2>>& va, const vector<NDPoint<2>>& vb) {
    const auto& a = va.front();
    const auto& b = vb.front();
    
    for (int i = 0; i < dim; ++i)
      boxes_[i].reset();
    
    if (dim == 1) {
      boxes_[1].reset(new DrawnBox(0, a[1], 1, b[1], "green"));      
      gogui::snapshot("PS Dim=1/2");
    } else if (dim == 0) {
      boxes_[0].reset(new DrawnBox(a[0], boxes_[1]->y1, b[0], boxes_[1]->y2, "green"));
      gogui::snapshot("PS Dim=0/2");
    }
  }
  
  void onLastDimFound(const vector<NDPoint<2>>& found) {
    for (const auto& p : found) {
      auto it = find(orig_.begin(), orig_.end(), gogui::Point{p[0], p[1]});
      if (it != orig_.end()) {
        it->setColor("orange");
      }
    }
    gogui::snapshot("Found a point");
  }
  
  gogui::vector<gogui::Point>& orig_;
  std::array<std::unique_ptr<DrawnBox>, 2> boxes_;
};

void demo() {  
  using O = ORT<2, vector<NDPoint<2>>, DimCmpInVec<2>, SumOfPointsTrans<2>>;
  
  gogui::vector<gogui::Point> points;
  vector<vector<NDPoint<2>>> pts;
  for (int i = 0; i < 350; ++i) {
    double x = static_cast<double>(rand()) / RAND_MAX;    
    double y = static_cast<double>(rand()) / RAND_MAX;
    pts.push_back({{x, y}});
    points.push_back({x, y});
  } 
  
  gogui::snapshot("Cloud");
  
  vector<NDPoint<2>> a = {{0.17, 0.26}};
  vector<NDPoint<2>> b = {{0.63, 0.87}};
  
  {
    gogui::Line l1{gogui::Point{a[0][0], a[0][1]}, gogui::Point{a[0][0], b[0][1]}}; l1.setColor("blue");
    gogui::Line l2{gogui::Point{a[0][0], b[0][1]}, gogui::Point{b[0][0], b[0][1]}}; l2.setColor("blue");
    gogui::Line l3{gogui::Point{b[0][0], b[0][1]}, gogui::Point{b[0][0], a[0][1]}}; l3.setColor("blue");
    gogui::Line l4{gogui::Point{b[0][0], a[0][1]}, gogui::Point{a[0][0], a[0][1]}}; l4.setColor("blue");
    gogui::ActiveLine al1{l1}, al2{l2}, al3{l3}, al4{l4};
    gogui::snapshot("Request");
    
    O tree(pts, SumOfPointsMix<2>{});
    bool any;
    {
      GoGuiVisualizer<2> viz{points};
      tree.query(a, b, any, viz);
    }
  }
  
  gogui::snapshot("Done");
}

int main(int argc, char** argv) {
  
  demo();
  std::ofstream("demo.json") << gogui::getJSON();
  gogui::reset();
  
  testConstructionForSumOfPoints<2>({10, 33, 100, 333, 1000, 3333, 10000, 33333, 100000}); 
  testConstructionForSumOfPoints<3>({10, 33, 100, 333, 1000, 3333, 10000}); 
  testConstructionForSumOfPoints<4>({10, 33, 100, 333, 1000, 3333});  
  testConstructionForSumOfPoints<5>({10, 33, 100, 333});   
  
  testConstructionForMax<2>({10, 100, 1000, 10000, 100000, 1000000});
  testConstructionForMax<3>({10, 100, 1000, 10000, 100000});
  testConstructionForMax<4>({10, 100, 1000, 10000});
  testConstructionForMax<5>({10, 100, 1000});
}
