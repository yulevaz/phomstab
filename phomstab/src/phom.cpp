#include<Rcpp.h>
#include <CGAL/Epick_d.h>
//[[Rcpp:depends(CGAL/Epick_d.h)]]
#include <gudhi/Simplex_tree.h>
//[[Rcpp:depends(gudhi/Simplex_tree.h)]]
#include <gudhi/Euclidean_witness_complex.h>
//[[Rcpp:depends(gudhi/Euclidean_witness_complex.h)]]
#include <gudhi/Euclidean_strong_witness_complex.h>
//[[Rcpp:depends(gudhi/Euclidean_strong_witness_complex.h)]]
#include <gudhi/pick_n_random_points.h>
//[[Rcpp:depends(gudhi/pick_n_random_points.h)]]
#include <gudhi/Persistent_cohomology.h>
//[[Rcpp:depends(gudhi/Persistent_cohomology)]]
#include <vector>

using namespace std;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef typename K::Point_d Point_d;
typedef typename Gudhi::witness_complex::Euclidean_strong_witness_complex<K> Strong_witness_complex;
typedef typename Gudhi::witness_complex::Euclidean_witness_complex<K> Weak_witness_complex;
typedef typename Gudhi::Simplex_tree<> Simplex_tree; 
typedef std::vector< Point_d > Point_vector;

//Persistence diagrams to Rcpp::Matrix (extracted from TDA package of CRAN: https://cran.r-project.org/web/packages/TDA/index.html)
/** \brief This function converts the StlMatrix objects to Rcpp::RcppMatrix
 * @param[out] 	Rcpp::RcppMatrix	A matrix containing the persistence diagram
 * @param[in]	sltMatrices		A std matrix containing the persistence diagram
 * @param[in]	includeIndex		Boolean value which includes elements indices on data if true
 * @param[in]	colNum			Number of columns in sltMatrices
 */
template<typename RcppMatrix, typename StlMatrix>
RcppMatrix concatStlToRcpp(const std::vector< StlMatrix >& stlMatrices,
		bool includeIndex, unsigned colNum) {
	unsigned rowNum = 0;

	typename std::vector< StlMatrix >::const_iterator vecItr;
	for (vecItr = stlMatrices.begin(); vecItr != stlMatrices.end(); ++vecItr) {
		rowNum += vecItr->size();
	}
	RcppMatrix rcppMatrix(rowNum, colNum);

	unsigned vecIdx, rowIdx, colIdx;
	for (vecIdx = 0, rowIdx = 0; vecIdx < stlMatrices.size(); ++vecIdx) {
		typename StlMatrix::const_iterator matItr;
		for (matItr = stlMatrices[vecIdx].begin();
				matItr != stlMatrices[vecIdx].end(); ++matItr, ++rowIdx) {
			if (includeIndex) {
				rcppMatrix[rowIdx] = vecIdx;
				for (colIdx = 0; colIdx < colNum - 1; ++colIdx) {
					rcppMatrix[rowIdx + (colIdx + 1) * rowNum] = (*matItr)[colIdx];
				}
			}
			else {
				for (colIdx = 0; colIdx < colNum; ++colIdx) {
					rcppMatrix[rowIdx + colIdx * rowNum] = (*matItr)[colIdx];
				}
			}
		}
	}

	return rcppMatrix;
}

//Convert data input to std::vector
/** \brief This function converts the data input provided by R to a std::vector object
 * @param[out] 	std::vector<Point_t>	A vector which contains the input data points
 * @param[in] 	X			The input data points provided by R
 */
Point_vector DataToVector(const Rcpp::NumericMatrix & X) {

	int nr = X.nrow();
	int nc = X.ncol();

	Point_vector data;

	for (int i = 0; i < nr; i++) {

		Point_d point;

		for (int j = 0; j < nc; j++) {

			double el = X(i,j);
			point.push_back(el);

		}

		data.push_back(point);

	}

	return data;

}

//Weak Witness from GUDHI package
/** \brief This function builds, on top of an input data points, a filtration employing lazy witness complex from GUDHI package. Note that the landmarks is chosen by random subsampling and the distanceadopted is the euclidean.
 * @param[out] Gudhi::Simplex_tree<>	A GUDHI simplex tree
 * 
 * @param[in] X			Input data points matrix
 * @param[in] nLandmarks	Number of landmarks adopted in lazy witness which are chosen by random subsampling
 * @param[in] maxDimension	Maximum dimesion on the constructed complex
 * @param[in] relaxed		If true the lazy witness is relaxed as shown in http://gudhi.gforge.inria.fr/doc/latest/group__witness__complex.html and alpha applies. 
 * @param[in] alpha		Parameter alpha of relaxed lazy witness
 */
Gudhi::Simplex_tree<> strongWitTree(const Rcpp::NumericMatrix & X
				, const int	nLandmarks
				, const int	maxDimension
				, const bool	relaxed
				, const double	alpha) {

	using Filtration_value = Simplex_tree::Filtration_value;

	Simplex_tree simplex_tree;
	Point_vector data, landmarks;
	//converting R input to std::vector<Point_d>
	data = DataToVector(X);

	double alp;

	if (!relaxed)
		alp = std::numeric_limits<Filtration_value>::infinity();
	else 
		alp = alpha;

	//creating simplex tree
	Gudhi::subsampling::pick_n_random_points(data, nLandmarks, std::back_inserter(landmarks));
	Strong_witness_complex witness_complex(landmarks, data);
	witness_complex.create_complex(simplex_tree, alp, maxDimension);
	simplex_tree.initialize_filtration();

	return simplex_tree;
}

//Strong Witness from GUDHI package
/** \brief This function builds, on top of an input data points, a filtration employing lazy witness complex from GUDHI package. Note that the landmarks is chosen by random subsampling and the distanceadopted is the euclidean
 * @param[out] Gudhi::Simplex_tree<>	A GUDHI simplex tree
 * 
 * @param[in] X			Input data points matrix
 * @param[in] nLandmarks	Number of landmarks adopted in lazy witness which are chosen by random subsampling
 * @param[in] maxDimension	Maximum dimesion on the constructed complex
 * @param[in] relaxed		If true the lazy witness is relaxed as shown in http://gudhi.gforge.inria.fr/doc/latest/group__witness__complex.html and alpha applies. 
 * @param[in] alpha		Parameter alpha of relaxed lazy witness
 */
Gudhi::Simplex_tree<> weakWitTree(const Rcpp::NumericMatrix & X
				, const int	nLandmarks
				, const int	maxDimension
				, const bool	relaxed
				, const double	alpha) {

	using Filtration_value = Simplex_tree::Filtration_value;

	Simplex_tree simplex_tree;
	Point_vector data, landmarks;
	//converting R input to std::vector<Point_d>
	data = DataToVector(X);

	double alp;

	if (!relaxed)
		alp = std::numeric_limits<Filtration_value>::infinity();
	else 
		alp = alpha;

	//creating simplex tree
	Gudhi::subsampling::pick_n_random_points(data, nLandmarks, std::back_inserter(landmarks));
	Weak_witness_complex witness_complex(landmarks, data);
	witness_complex.create_complex(simplex_tree, alp, maxDimension);
	simplex_tree.initialize_filtration();

	return simplex_tree;
}

// Persistence diagram builder
/** \brief This function produces and returns the persistence diagram
 * @param[out]	Rcpp::List	A R list type persistence diagram
 * @param[in] 	radiusInit	The initial radius of the filtration
 * @param[in] 	maxDimension	Maximum dimesion on the constructed complex
 * @param[in] 	simplex_tree	The simplex tree built by strong or weak witness complex
 * @param[out]	Rcpp::List	A R list type persistence diagram
 */
Rcpp::List defineIntervals(int radiusInit, int maxDimension, Simplex_tree simplex_tree) {
    	
	using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
      	using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
	//simplex_tree.print_hasse(cout);
	std::vector<std::vector<std::vector<double>>> persDgm;

	//persistence diagram
	Persistent_cohomology pcoh(simplex_tree);
	pcoh.init_coefficients(2);
	pcoh.compute_persistent_cohomology(radiusInit);

	pcoh.output_diagram();
        //extracted from TDA R package (https://cran.r-project.org/web/packages/TDA/index.html)	gudhiUtils.h / FiltrationDiagGudhi 
	std::vector<double> dgmPoint(2);
  	std::vector<std::vector<double>> dgm = pcoh.memory_output_diagram<std::vector<std::vector<double>>>();

	//printDiag(dgm);	

	persDgm.resize(maxDimension);
	for (unsigned rowIdx = 0; rowIdx < dgm.size(); ++rowIdx) {
		dgmPoint[0] = dgm[rowIdx][2];
		dgmPoint[1] = dgm[rowIdx][3];
		persDgm[dgm[rowIdx][1]].push_back(dgmPoint);
	}
	////////////////////////////////////////////////////////

	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3));
}

//' This function builds, on top of an input data points, a filtration employing weak or strong witness complex from GUDHI package. Note that the landmarks is chosen by random subsampling and the distance adopted is the euclidean.
//'
//' param X		Input data points matrix
//' param nLandmarks	Number of landmarks adopted in lazy witness which are chosen by random subsampling
//' param maxDimension	Maximum dimesion on the constructed complex
//' param isWeak	Defines if the program should employ weak witness or strong witness
//' param relaxed	If true the lazy witness is relaxed as shown in http://gudhi.gforge.inria.fr/doc/latest/group__witness__complex.html and alpha applies. 
//' param alpha		Parameter alpha of relaxed lazy witness
//' return		A persistence diagram
//' @export
//[[Rcpp::export]]
Rcpp::List witDiag(const Rcpp::NumericMatrix & X
				, const int	nLandmarks
				, const int	maxDimension
				, const bool	isWeak
				, const int	relaxed
				, const double	alpha) {

	if (isWeak) {
		Simplex_tree simplex_tree = weakWitTree(X,nLandmarks,maxDimension,relaxed,alpha);
		return defineIntervals(-1, maxDimension, simplex_tree);
	}
	else {
		Simplex_tree simplex_tree = strongWitTree(X,nLandmarks,maxDimension,relaxed,alpha);
		return defineIntervals(-1, maxDimension, simplex_tree);
	}

}

