namespace Corr {
class TimeSliceFunc : public SetFunc {
public:
  TimeSliceFunc(int dir) : dir_decay(dir) {}

  int operator()(const multi1d<int> &coordinate) const;

  int numSubsets() const;

private:
  TimeSliceFunc() {} // hide default constructor

  int dir_decay;
};

int TimeSliceFunc::operator()(const multi1d<int> &coordinate) const {
  if ((dir_decay < 0) || (dir_decay >= Nd)) {
    return 0;
  } else {
    return coordinate[dir_decay];
  }
}

int TimeSliceFunc::numSubsets() const {
  if ((dir_decay < 0) || (dir_decay >= Nd)) {
    return 1;
  } else {
    return Layout::lattSize()[dir_decay];
  }
}
}
