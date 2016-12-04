template<size_t Dim>
using NDPoint = std::array<double, Dim>;

template<size_t Dim>
struct DimCmpSingle {
  template<size_t IthDim>
  static bool precedes(
      const NDPoint<Dim>& p1,
      const NDPoint<Dim>& p2) {
    return p1[IthDim] < p2[IthDim];
  }
};

template<size_t Dim>
struct MaxFn {
  NDPoint<Dim+1> operator()(
      const NDPoint<Dim+1>& p1,
      const NDPoint<Dim+1>& p2) const {
    return p1[Dim] > p2[Dim] ? p1 : p2;
  }
};

auto createMax3DTree(const std::vector<NDPoint<4>>& data) {
  return ORT<3, NDPoint<4>, DimCmpSingle<3>> tree(data, MaxFn<3>{});
}
