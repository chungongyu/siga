#ifndef runner_h_
#define runner_h_

#include <config.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <getopt.h>

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>

typedef boost::property_tree::ptree Properties;
typedef std::vector<std::string> Arguments;

class Runner {
 public:
  const std::string& options(const option** longopts = nullptr) const {
    if (longopts != nullptr) {
      *longopts = _longopts;
    }
    return _shortopts;
  }
  std::string transform(char key, int* flag = nullptr) const {
    if (_longopts != nullptr) {
      for (size_t i = 0; _longopts[i].name != nullptr; ++i) {
        if (key == _longopts[i].val) {
          if (flag != nullptr) {
            *flag = _longopts[i].has_arg;
          }
          return _longopts[i].name;
        }
      }
    }
    return std::string(1, key);
  }
  virtual int run(const Properties& options, const Arguments& arguments) = 0;

 protected:
  Runner(const std::string& shortopts = "", const option* longopts = nullptr)
      : _shortopts(shortopts), _longopts(longopts) {
  }
  std::string _shortopts;
  const option* _longopts;
};

typedef Runner* RunnerPtr;

class RunnerManager {
 public:
  RunnerManager() {
  }
  virtual ~RunnerManager() {
  }

  static RunnerManager* get() {
    static RunnerManager mgr;
    return &mgr;
  }

  RunnerPtr create(const std::string& name) const {
    auto i = _runners.find(name);
    if (i != _runners.end()) {
      return std::get<0>(i->second);
    }
    return RunnerPtr();
  }
  bool install(const std::string& name, RunnerPtr runner, const std::string& introduction, int rank = 0) {
    if (_runners.find(name) != _runners.end()) {
      return false;
    }
    _runners[name] = std::make_tuple(runner, introduction, rank);
    return true;
  }
  bool uninstall(const std::string& name) {
    auto i = _runners.find(name);
    if (i != _runners.end()) {
      _runners.erase(i);
      return true;
    }
    return false;
  }

  int help(int ret = 0) const {
    std::cout << boost::format("%s version %s, report bugs to [%s]")
        % PACKAGE_NAME % PACKAGE_VERSION % PACKAGE_BUGREPORT << std::endl;
    std::cout << std::endl;
    std::cout << boost::format("Usage: %s <command> [<args>]") % PACKAGE_NAME << std::endl;
    std::cout << std::endl;
    std::cout << boost::format("The most commonly used %s commands are:") % PACKAGE_NAME << std::endl;

    {
      size_t max_name_length = 0;
      typedef std::vector<std::tuple<std::string, int> > _RunnerList_;
      _RunnerList_ runners;
      {
        for (auto i = _runners.begin(); i != _runners.end(); ++i) {
          if (std::get<2>(i->second) != kUnknown) {
            max_name_length = std::max(max_name_length, i->first.length());
            runners.push_back(std::make_tuple(i->first, std::get<2>(i->second)));
          }
        }
      }
      // sort by runner's rank
      std::sort(runners.begin(), runners.end(),
          [&](const std::tuple<std::string, int>& l, const std::tuple<std::string, int>& r) {
            auto x = std::get<1>(l), y = std::get<1>(r);
            if (x != y) {
              return x < y;
            }
            return std::get<0>(l) < std::get<0>(r);
          });
      max_name_length += 2;
      for (auto i = runners.begin(); i != runners.end(); ++i) {
        if (std::get<1>(*i) != kUnknown) {
          std::string cmd = std::get<0>(*i);
          auto info = _runners.find(cmd);
          assert(info != _runners.end());
          cmd.resize(max_name_length, ' ');
          std::cout << boost::format("   %s%s") % cmd % std::get<1>(info->second) << std::endl;
        }
      }
    }

    std::cout << std::endl;
    std::cout << boost::format("See '%s <command> -h' to read about a specific subcommand.")
        % PACKAGE_NAME << std::endl;
    std::cout << boost::format("Further help: %s") % PACKAGE_URL << std::endl;

    return ret;
  }
  int help(int argc, char* argv[]) const {
    static std::string shortopt("vh");
    static option longopts[] = {
      {"version",   no_argument, nullptr, 'v'},
      {"help",    no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}
    };
    int opt = -1;
    while ((opt = getopt_long(argc, argv, shortopt.c_str(), longopts, nullptr)) != -1) {
      switch ((char)opt) {
      case 'v':
        std::cout << boost::format("%s version %s") % PACKAGE_NAME % PACKAGE_VERSION << std::endl;
        return 0;
      case 'h':
        return help(0);
      }
    }
    return help(1);
  }

 private:
  typedef std::tuple<RunnerPtr, std::string, int> RunnerInfo;
  typedef std::map<std::string, RunnerInfo> RunnerList;
  RunnerList _runners;
};

#define RUNNER_INSTALL(name, runner, introduction, rank) \
  RunnerManager::get()->install(name, runner, introduction, rank)
#define RUNNER_UNINSTALL(name) \
  RunnerManager::get()->uninstall(name)

#endif  // runner_h_
