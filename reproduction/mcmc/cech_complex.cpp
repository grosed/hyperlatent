#include <Rcpp.h>
using namespace Rcpp;
#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>

// [[Rcpp::export]]

List calc_cech_R2(float rmax, NumericMatrix xs, int nrow, int ncol, int K) {
  // std::cout << nrow << " ";
  // std::cout << ncol << " ";
  // std::cout << std::endl;
  // Type definitions
  using Point_cloud = std::vector<std::array<double, 2>>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;
  Point_cloud points;

  // ----------------------------------------------------------------------------
  // Populate points with xs
  // ----------------------------------------------------------------------------
  std::array<double, 2> xtmp; 
  for (int i=0; i < nrow; i++){
    for (int j=0; j < ncol; j++ ){
      xtmp[j] = xs[i*ncol + j];
      //std::cout << xtmp[j] << " ";
    }
    //std::cout << std::endl;
    points.push_back(xtmp); // adds element to the end of a vector
  }

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from points
  // ----------------------------------------------------------------------------
  Filtration_value max_radius = rmax;
  Cech_complex cech_complex_from_points(points, max_radius);
  Simplex_tree stree;
  cech_complex_from_points.create_complex(stree, K); // calculate complex up to order K interaction 
  // ----------------------------------------------------------------------------
  // Extract information about the one skeleton Cech complex
  // ----------------------------------------------------------------------------
  int n = stree.num_simplices(); // number of simplices 
  List out_nodes(n); // list of length n
  List out_filt(n); // list of length n
  int i = 0;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    out_nodes[i] = stree.simplex_vertex_range(f_simplex);
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      out_filt[i] = stree.filtration(f_simplex);
    }
    i ++ ;
  }

  return List::create(_["nodes"]=out_nodes, _["filts"]=out_filt);
}

// [[Rcpp::export]]

List calc_cech_R3(float rmax, NumericMatrix xs, int nrow, int ncol, int K) {
  // std::cout << nrow << " ";
  // std::cout << ncol << " ";
  // std::cout << std::endl;
  // Type definitions
  using Point_cloud = std::vector<std::array<double, 3>>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;
  Point_cloud points;

  // ----------------------------------------------------------------------------
  // Populate points with xs
  // ----------------------------------------------------------------------------
  std::array<double, 3> xtmp; 
  for (int i=0; i < nrow; i++){
    for (int j=0; j < ncol; j++ ){
      xtmp[j] = xs[i*ncol + j];
      //std::cout << xtmp[j] << " ";
    }
    //std::cout << std::endl;
    points.push_back(xtmp); // adds element to the end of a vector
  }

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from points
  // ----------------------------------------------------------------------------
  Filtration_value max_radius = rmax;
  Cech_complex cech_complex_from_points(points, max_radius);
  Simplex_tree stree;
  cech_complex_from_points.create_complex(stree, K);
  // ----------------------------------------------------------------------------
  // Extract information about the one skeleton Cech complex
  // ----------------------------------------------------------------------------
  int n = stree.num_simplices(); // number of simplices 
  List out_nodes(n); // list of length n
  List out_filt(n); // list of length n
  int i = 0;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    out_nodes[i] = stree.simplex_vertex_range(f_simplex);
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      out_filt[i] = stree.filtration(f_simplex);
    }
    i ++ ;
  }

  return List::create(_["nodes"]=out_nodes, _["filts"]=out_filt);
}
