#include <gtest/gtest.h>

#include <iostream>

int main(int argc, char **argv) {
  srand(time(0));
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "RUNNING TESTS ..." << std::endl;
  int ret{RUN_ALL_TESTS()};
  if (!ret)
    std::cout << "<<<SUCCESS>>>" << std::endl;
  else
    std::cout << "FAILED" << std::endl;
  return 0;
}