#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <fstream>
#include <iterator>
#include "math.h"

#define MAXBUFSIZE  ((int) 1e6)

using Eigen::MatrixXd;

class Point2D
{
	friend class FElement;

	float x, y; //PunktKoordinaten
	int u, v; // Nummer der Freiheitsgrade

public:
	Point2D()
	{
		x = y = u = v = 0;
	}
	Point2D(float X, float Y) : x(X), y(Y)
	{
		u = v = 0;
	};
	void setX(float X)
	{
		x = X;
	}
	void setY(float Y)
	{
		y = Y;
	}
	void setPoint(float X, float Y)
	{
		x = X;
		y = Y;
	}
	float getX() const
	{
		return x;
	}
	float getY() const
	{
		return y;
	}
};

class FElement
{
	friend class Struktur;
	Point2D SP, EP; 
	float E, A, L;
	Eigen::MatrixXf KMAT;
	std::vector<size_t> edofVec;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		//Konstruktoren
	FElement()
	{
		SP.setPoint(0, 0);
		EP.setPoint(0, 0);
		E = A = L = 0;
	};
	FElement(float XSP, float YSP, float XEP, float YEP, float e, float a) : E(e), A(a)
	{
		SP.setPoint(XSP, YSP);
		EP.setPoint(XEP, YEP);
		L = calcL(SP, EP);
		KMAT = calcKmat(XSP, YSP, XEP, YEP, e, a, L);
	};
	FElement(Point2D* sp, Point2D* ep, float e, float a) : E(e), A(a)
	{
		SP.x = sp->x;
		SP.y = sp->y;
		EP.x = ep->x;
		EP.y = ep->y;
		L = calcL(SP, EP);
		KMAT = calcKmat(SP, EP, e, a, L);
	};
	FElement(Point2D& sp, Point2D& ep, float e, float a) : E(e), A(a)
	{
		SP = sp;
		EP = ep;
		L = calcL(SP, EP);
		KMAT = calcKmat(SP, EP, e, a, L);
	};

	float calcL(Point2D& SP, Point2D& EP)
	{
		//Berechnet die Elementlänge ausgehend vom Start- und Endpunkt
		L = sqrt(pow(EP.getX() - SP.getX(), 2) + pow(EP.getY() - SP.getY(), 2));
		return L;
	}

	Eigen::Matrix4f calcKmat(Point2D& SP, Point2D& EP, float e, float a, float L)
	{
		Eigen::MatrixXf T(2, 4), TTrans(4, 2);
		Eigen::Matrix2f K;
		T << EP.getX() - SP.getX(), EP.getY() - SP.getY(), 0, 0,
			0, 0, EP.getX() - SP.getX(), EP.getY() - SP.getY();
		T = 1 / L * T;

		TTrans = T.transpose();

		K << 1, -1, -1, 1;
		K = e * a / L * K;

		KMAT = TTrans * K * T;
		return KMAT;
	}

	Eigen::Matrix4f calcKmat(float XSP, float YSP, float XEP, float YEP, float e, float a, float L)
	{
		Eigen::MatrixXf T(2, 4), TTrans(4, 2);
		Eigen::Matrix2f K;
		T << XEP - XSP, YEP - YSP, 0, 0,
			0, 0, XEP - XSP, YEP - YSP;
		T = 1 / L * T;

		TTrans = T.transpose();

		K << 1, -1, -1, 1;
		K = e * a / L * K;

		KMAT = TTrans * K * T;
		return KMAT;
	}

	void printKmat() const
	{
		std::cout << KMAT << std::endl;
	}
	Eigen::MatrixXf retKmat() const
	{
		return KMAT;
	}

	void setSPU(int U)
	{
		SP.u = U;
	}
	void setSPV(int V)
	{
		SP.v = V;
	}
	void setEPU(int U)
	{
		EP.u = U;
	}
	void setEPV(int V)
	{
		EP.v = V;
	}

	void setSPPoint(Point2D& sp)
	{
		SP = sp;
	}
	void setEPPoint(Point2D& ep)
	{
		EP = ep;
	}

	Point2D getSP() const
	{
		return SP;
	}
	Point2D getEP() const
	{
		return EP;
	}

	int getSPU() const
	{
		return SP.u;
	}
	int getSPV() const
	{
		return SP.v;
	}
	int getEPU() const
	{
		return EP.u;
	}
	int getEPV() const
	{
		return EP.v;
	}

	friend std::ostream& operator <<(std::ostream& os, const FElement& FEEL)
	{
		os << "Querschnittsfläche: " << FEEL.A << std::endl;
		os << "Elastizitätsmodul: " << FEEL.E << std::endl;
		os << "Elementlänge: " << FEEL.L << std::endl;
		os << "x-Koordinate - Startpunkt: " << FEEL.SP.getX() << std::endl;
		os << "y-Koordinate - Startpunkt: " << FEEL.SP.getY() << std::endl;
		os << "x-Koordinate - Endpunkt: " << FEEL.EP.getX() << std::endl;
		os << "y-Koordinate - Endpunkt: " << FEEL.EP.getY() << std::endl;
		os << "Nummer des Freiheitsgrads u - Startpunkt: " << FEEL.getSPU() << std::endl;
		os << "Nummer des Freiheitsgrads v - Startpunkt: " << FEEL.getSPV() << std::endl;
		os << "Nummer des Freiheitsgrads u - Endpunkt: " << FEEL.getEPU() << std::endl;
		os << "Nummer des Freiheitsgrads v - Endpunkt: " << FEEL.getEPV() << std::endl;
		os << std::endl;

		return os;
	}
};

class Struktur
{
	//Anpassen an fixed KMAT
	std::vector<FElement> StrukturVektor;
	Eigen::MatrixXf strctMat;
	Eigen::MatrixXf topologyMat;
	Eigen::VectorXf forceVec;
	std::vector<int> boundVec;

public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	void initSystem(std::string filepath)
	{
		std::ifstream dataStream;
		dataStream.open(filepath, std::ios_base::in);
		std::vector<std::vector<float>> dataVec;
		if (dataStream.is_open())
		{
			std::string lineBuffer;
			while (getline(dataStream, lineBuffer))
			{
				std::istringstream is(lineBuffer);
				dataVec.push_back(std::vector<float>(std::istream_iterator<float>(is),
				                                     std::istream_iterator<float>()));
			}
		}
		else
		{
			std::cerr << "Couldn't open the file" << std::endl;
		}

		for (int i = 0; i < dataVec.size(); i++)
		{
			FElement* felement = new FElement(dataVec[i][0],
			                                  dataVec[i][1],
			                                  dataVec[i][2],
			                                  dataVec[i][3],
			                                  dataVec[i][4],
			                                  dataVec[i][5]);
			StrukturVektor.push_back(*felement);
		}
	}

	void setNrOfEDOF()
	{
		//Erstes Element ist das Grundelement
		StrukturVektor[0].setSPU(1);
		StrukturVektor[0].setSPV(2);
		StrukturVektor[0].setEPU(3);
		StrukturVektor[0].setEPV(4);

		for (int i = 1; i < StrukturVektor.size(); ++i)
		{
			StrukturVektor[i].setSPPoint(StrukturVektor[i - 1].getEP());
			StrukturVektor[i].setEPU(StrukturVektor[i].getSPU() + 2);
			StrukturVektor[i].setEPV(StrukturVektor[i].getSPV() + 2);
		}
	}

	void compileTopology()
	{
		topologyMat.resize(StrukturVektor.size(), 5);

		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			topologyMat.col(0)(i, 0) = i;
		}		
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			topologyMat.col(1)(i, 0) = StrukturVektor[i].getSPU();
		}
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			topologyMat.col(2)(i, 0) = StrukturVektor[i].getSPV();
		}
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			topologyMat.col(3)(i, 0) = StrukturVektor[i].getEPU();
		}
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			topologyMat.col(4)(i, 0) = StrukturVektor[i].getEPV();
		}
	}

	void readTopologyFromFile(std::string filepath)
	{
		
		std::ifstream fin(filepath);
		topologyMat.resize(3, 5);


		if (fin.is_open())
		{
			for (int row = 0; row < 3; row++)
				for (int col = 0; col < 5; col++)
				{
					float item = 0.0;
					fin >> item;
					topologyMat(row, col) = item;
				}
			fin.close();
		}
	
	}


	void compileStructureMatrix()
	{
		//Testen, ob Topology existiert

		//Anlegen der Topology-Untermatrix
		Eigen::MatrixXf TUntermatrix;
		TUntermatrix = topologyMat.block(0, 1, topologyMat.rows(), 4);
		strctMat = Eigen::MatrixXf::Zero((StrukturVektor.size() + 1) * 2, (StrukturVektor.size() + 1) * 2);

		//Äußere Iteration über die Topologie-Untermatrix
		for (int i = 0; i < TUntermatrix.rows(); ++i)
		{
			//TUntermatrix.col(i).size()
			//Innere Iteration über die Elemte der Topologie-Untermatrix
			for (int j = 0; j < 4; ++j)
			{
				for (int k = 0; k < 4; ++k)
				{
					strctMat(TUntermatrix(i, j) - 1, TUntermatrix(i, k) - 1) += StrukturVektor[i].KMAT(j, k);
					std::cout << strctMat << std::endl << std::endl;
				}
			}
		}
	}

	void setBoundaries(std::string filepath)
	{
		std::ifstream fin(filepath);

		int num;

		if (fin.is_open())
		{
			//Testen ob der Kraftvektor vollständig gefüllt ist ! Sonst 0 setzen -> throw

			while (fin >> num)
				boundVec.push_back(num);

			fin.close();
		}
		//Setze Spalte und Zeile gleich 0 	
		//und das Diagonalelement zu 1
		for (int i = 0; i < boundVec.size(); ++i)
		{
			strctMat.col(boundVec[i] - 1).setZero();
			strctMat.row(boundVec[i] - 1).setZero();
			strctMat(boundVec[i] - 1, boundVec[i] - 1) = 1;
			forceVec(boundVec[i] - 1) = 0;
		}
	}

	void setForceVec(std::string filepath)
	{
		forceVec.resize((StrukturVektor.size() + 1) * 2);
		std::ifstream fin(filepath);

		if (fin.is_open())
		{
			//Testen ob der Kraftvektor vollständig gefüllt ist ! Sonst 0 setzen -> throw

			for (int i = 0; i < (StrukturVektor.size() + 1) * 2; i++)
			
				fin >> forceVec[i];

			fin.close();
		}

		for(int i = 0 ; i < forceVec.size();++i)
		{
			forceVec(i) *= -1;
		}
	}

	void solve()
	{
		Eigen::VectorXf x = strctMat.householderQr().solve(forceVec);
		std::cout << x;
		std::cout << std::endl;
	}

	//Hilfsfunktionen
	void printMatrices()
	{
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			std::cout << "Elementmatrix: " << i+1 << std::endl;
			StrukturVektor[i].printKmat();
			std::cout << std::endl;
		}
	}

	void printMatricesInFile(std::ofstream& fout)
	{
		fout.open("Elementmatrizen.txt", std::ios::out);
		if (!fout.is_open())
		{
			std::cerr << "open error\n";
		}
		for (int k = 0; k < StrukturVektor.size(); ++k)
		{
			fout << "Elementmatrix: " << k+1 << std::endl;
			fout << StrukturVektor[k].retKmat();
			fout << std::endl;
		}
		fout.close();
	}

	void printCompiledMatrix()
	{
		std::cout << std::endl;
		std::cout << "Gesamtmatrix des Systems: " << std::endl;
		std::cout << strctMat;
		std::cout << std::endl;
	}

	void printForceVec()
	{
		std::cout << std::endl;
		std::cout << "Kraftvektor: " << std::endl;
		std::cout << forceVec;
		std::cout << std::endl;
	}

	void printAllElements()
	{
		for (int i = 0; i < StrukturVektor.size(); ++i)
		{
			std::cout << "Element: " << i+1 << std::endl;
			std::cout << StrukturVektor[i] << std::endl;
		}
	}

	void printTopology()
	{
		std::cout << std::endl;
		std::cout << "Topologiematrix des Systems: " << std::endl;
		std::cout << topologyMat;
		std::cout << std::endl;
	}
};

int main(int argc, char** argv)
{
	std::ofstream fileout;
	Struktur * Stabwerk = new Struktur;
	Stabwerk->initSystem(argv[1]);
	Stabwerk->setForceVec("D:\\Programme\\Git\\Repos\\FEM2D\\FEM2D\\kraftvektor.txt");
	Stabwerk->setNrOfEDOF();
	Stabwerk->printMatrices();
	Stabwerk->printAllElements();
	//Stabwerk->compileTopology();
	Stabwerk->readTopologyFromFile("D:\\Programme\\Git\\Repos\\FEM2D\\FEM2D\\topology.txt");
	Stabwerk->printTopology();
	Stabwerk->compileStructureMatrix();
	Stabwerk->printCompiledMatrix();
	Stabwerk->setBoundaries("D:\\Programme\\Git\\Repos\\FEM2D\\FEM2D\\boundaries.txt");
	Stabwerk->printForceVec();
	Stabwerk->printCompiledMatrix();
	Stabwerk->solve();
}