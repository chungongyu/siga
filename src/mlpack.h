#ifndef mlpack_h_
#define mlpack_h_

#include <memory>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/serialization/serialization.hpp>

#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>

//
// This is the class that we will serialize.  It is a pretty simple wrapper
// around mlpack Model.
//

template <typename Model>
class AIModel {
 public:
  // The tree itself, left public for direct access by this program.
  Model model;
  mlpack::data::DatasetInfo info;

  // Serialize the model.
  template <typename Archive>
  void serialize(Archive& ar, const unsigned int /* version */) {
    ar & BOOST_SERIALIZATION_NVP(model);
    ar & BOOST_SERIALIZATION_NVP(info);
  }

  bool load(const std::string& filename) {
    return mlpack::data::Load(filename, "model", *this, false);
  }

  template <typename VecType>
  size_t Classify(const VecType& vec) const {
    return model.Classify(vec);
  }
};

template <typename Model>
class BaggingModel {
 public:
  bool load(const std::vector<std::string>& filelist) {
    for (const auto& filename : filelist) {
      AIModelPtr model = AIModelPtr(new AIModel<Model>());
      if (!model->load(filename)) {
        return false;
      }
      _models.push_back(model);
    }
    return true;
  }
  bool load(const std::string& filelist) {
    std::vector<std::string> vec;
    boost::algorithm::split(vec, filelist, boost::algorithm::is_any_of(","));
    return load(vec);
  }

  template <typename VecType>
  size_t Classify(const VecType& vec) const {
    size_t l = 0;
    for (auto m : _models) {
      l += m->Classify(vec);
    }
    return l;
  }

  size_t Size() const {
    return _models.size();
  }

 private:
  typedef std::shared_ptr<AIModel<Model> > AIModelPtr;
  std::vector<AIModelPtr > _models;
};

#endif  // mlpack_h_
