#ifndef LABELPAIR_HPP
#define LABELPAIR_HPP

#include "unap.hpp"

namespace UNAP
{
class labelPair
{
private:
  label first_;
  label second_;
  label faceI_;

public:
  labelPair(label i, label j);

  inline label first() const { return first_; }

  inline label second() const { return second_; }

  inline label faceI() const { return faceI_; }

  inline void first(const label i) { first_ = i; }

  inline void second(const label i) { second_ = i; }

  inline void faceI(const label i) { faceI_ = i; }

  inline bool operator==(const labelPair &a) const
  {
    return ((this->first() == a.first()) && (this->second() == a.second()));
  }
};

}  // namespace UNAP

#endif  //- LABELPAIR_HPP
