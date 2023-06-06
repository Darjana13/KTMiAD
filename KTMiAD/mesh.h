#pragma once

#include "includes.h"
#include "Gauss.h"
#include "cgm.h"
#include "Solver.h"


struct type_points // ����� �������� 
{ 
public:
   double x, y, z; 

   type_points GetPoint(type_points b)
   {
	   type_points p;
	   p.x = (x + b.x) / 2;
	   p.y = (y + b.y) / 2;
	   p.z = (z + b.z) / 2;
	   return p;
   }
};

class c_line //����� ����� ��������� 
{
public:
   vector <double> points; 
};

class type_coef_lines
{
public:
   //c_line& start, end;
   int start, end;
   int n;
   double k;

   int SplitSegment(vector<c_line>& base,vector<c_line> &result, int size_points); // ��������� �� ������������ �������� �����
   //int SplitSegment(vector<c_line>& base, vector<c_line>& result); // ��������� �� ������������ �������� �����
   //int SplitSegment_z(vector<c_line>& base, vector<c_line>& result); // ��������� �� ������������ �������� �����

};

class type_elem // ����� ��������
{
public:
   vector <int> node_loc; // ������������ ��������� � ���������� ��������� �����
   vector <int> face_loc; // ������������ ��������� � ���������� ��������� ������
   vector <int> edge_loc; // ������������ ��������� � ���������� ��������� �����
   vector <int> bf_loc; // ������������ ��������� � ���������� ��������� �������� �������

   int obl; // ����� ������� �� �� ����� ���������
};

class type_area // ������� 
{
public:
   int material; //������ ���������
   int x_start; //������ �� �
   int x_end; // ����� �� �
   int y_start;
   int y_end;
   int z_start;
   int z_end;
};

class type_material
{
public:
	double lambda = 1;
	double gamma = 1;
};

struct BorderCond
{
	int elem, loc_fc;
	double value;

	const bool operator<(const BorderCond& rv) const {
		return (elem < rv.elem);
	}
};

class mesh
{
public:
   vector <c_line> x_line, y_line, z_line; // ������������ ����� �� ��z ��� ����� , ������� �������
   vector <c_line> x_base_line, y_base_line, z_base_line; // ������������ ����� �� ��z ������� ��������, �� ������ �� �� �����

   vector <type_area> area; // ������ ��������
   //vector<type_coef_lines> coef; // ������������ ��������� ������ ����������
   vector<type_coef_lines> coef_x, coef_y, coef_z;

   map<int, int> glob2real; // glob2real[i] = k ������ ��� i ����� ������ ����� ������������� k ������ �������� �����
   map<int, int> face_glob2real; // glob2real[i] = k ������ ��� i ����� ������ ����� ������������� k ������ �������� �����
   map<int, int> edge_glob2real; // glob2real[i] = k ������ ��� i ����� ������ ����� ������������� k ������ �������� �����
   map<int, int> bf_glob2real; // glob2real[i] = k ������ ��� i ����� ������ ����� ������������� k ������ �������� �����

   vector<set<int>> face_to_elem; // ��������, ��������� � ������
   vector<set<int>> edge_to_elem; // ��������, ��������� � ������

   vector <type_elem> all_elems; //������ ���������
   vector <type_points> all_nodes;
   vector<type_material> all_material;
   set<pair<int, double>> S1;

   int z_size = 0, y_size = 0, x_size = 0;
   int eg_all, fc_all, bf_all;
   int num_el, num_nod; // ���������� ��������� � ���������� �����

   int FillArea2(); // ������� ��������� ����� ������ ��������� � ����������
   int AddNodeToRealNodes(int glob_num, int i_node); // 
   
   int ReadMesh();
   int ReadMat();
   int InitMemory();
   int FillArea();
   int BuildS1();

  
   int GetEdge(int p1, int p2);
   int GetNeibEdge(int edge, int& p1, int& p2, set<int> &elems);
   int GetFace(int p1, int p2, int p3, int p4);
   int GetNeibFace(int face, set<int>& nodes, set<int>& elems);
   type_points GetCoordBF(int el, int bf_loc);
   type_points GetCentreFace(int el, int fc_loc);
   type_points GetCentreElem(int el);

   int WriteMesh(string path);
   int SplitSegment(const type_coef_lines& coef, c_line& lines);

};

class SLAU
{
public:
	vector<int> ia, ja;
	//vector<double> di, a, b;
	vector<double> di, au, al, b;

	int GeneratePortret(mesh& mesh);
	int AddLocal(vector<vector<double>>& A_loc, vector<double>& b_loc, int el_id, mesh& mesh);
	int SetS1_no_simmetry(mesh& mesh);
	int SetS1_test(mesh& mesh); // ���������� ������ �� �� �����, � �������� ������� �� includes 

	int ShowMatrix();   // ����������� ������� ������, ������� � ��������
	void ClearValues(); // ��������� 0 di, a, b
};

class FEM
{
	vector<vector<double>> A_loc;
	vector<double> b_loc;
	vector<double> q;

	vector<basis_func> basis; // �������� ���������� �������; ksi = (x - x_min) / (x_max - x_min))
	vector<basis_func> derivative_basis; // ����������� �������� ���������� �������; ksi = (x - x_min) / (x_max - x_min))

	Integrate_Gauss3Method3D gauss3method_3D;

	int GetG_Loc(int el_id, double B = 0);
	int GetG_Loc_test(int el_id);

	int Getb_Loc(int el_id);
	int Getb_Loc_test(int el_id);

public:
	mesh mesh;
	SLAU slau;

	int Init(string path);

	//int SolveTask(string path);                  // ���������� �����
	int SolveTask_test(string path);             // ���������� ������ �� �� �����, � �������� ������� �� includes 

	int GetP(vector<double>& P);                    // �������� P �� ���� �����
	int WriteP(string path);                    // �������� P �� ���� �����


	double GetPInPoint(Node3D point, int cur_el = -1); // �������� P � ����� (����� ������� �������, ���� �� ��������)
};


