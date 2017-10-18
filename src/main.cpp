#include "file_input.cpp"
#include "sim_input.cpp"

int main(int argc, char** argv)
{
  if (argc > 1) return file_input(argc, argv); // if there are arguments supplied - use file input
  return sim_input();                          // if no arguments are supplied use simulator input
}
