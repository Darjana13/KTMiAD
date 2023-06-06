#include "mesh.h"

int FEM::Init(string path)
{
	basis = { linear1, linear2 };
	derivative_basis = { dlinear1_dksi, dlinear2_dksi };
	gauss3method_3D.Init(basis, derivative_basis);

	mesh.ReadMesh();
	mesh.InitMemory();
	mesh.FillArea2();
	mesh.BuildS1();

	slau.GeneratePortret(mesh);
	q.resize(mesh.num_nod);
	b_loc.resize(8);
	A_loc.resize(8);
	for (int i = 0; i < 8; i++)
	{
		A_loc[i].resize(8);
	}
	return 0;
}

int FEM::SolveTask_test(string path)
{
	ofstream file;
	for (int i = 0; i < mesh.num_el; i++)
	{
		GetG_Loc_test(i);
		Getb_Loc_test(i); // Add M, if it exist
		slau.AddLocal(A_loc, b_loc, i, mesh);
		/*file.open("matrixA_" + toString(i));
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				file << A_loc[i][j] << "\t";
			}
			file << " | " << b_loc[i];
			file << '\n';
		}
		file.close();
		file.clear();*/
	}
	cout << "get SLAU" << endl;
	//slau.SetS1_test(mesh);

	//cgm_solver solver(slau.ia, slau.ja, slau.di, slau.a, slau.b);
	//solver.no_preconditioning(q);


	/*Solver check(slau.ia, slau.ja, slau.di, slau.au, slau.al, slau.b);
	vector<double> b_ist(mesh.num_nod);
	vector<double> q_ist(mesh.num_nod);
	for (int i = 0; i < mesh.num_nod; i++)
	{
		q_ist[i] = u_test(mesh.all_nodes[i].x, mesh.all_nodes[i].y, mesh.all_nodes[i].z);
	}
	check.A.Ax(q_ist, b_ist);
	for (int i = 0; i < mesh.num_nod; i++)
	{
		if (abs(slau.b[i] - b_ist[i]) > 1e-5)
			cout << " i " << i << " b[i] " << slau.b[i] << " Aq[i] " << b_ist[i] << endl;
	}*/

	slau.SetS1_no_simmetry(mesh);

	Solver solver(slau.ia, slau.ja, slau.di, slau.au, slau.al, slau.b);
	cout << "start solve " << endl;
	solver.CGM_LU();
	solver.getx0(q);



	/*Solver check(slau.ia, slau.ja, slau.di, slau.au, slau.al, slau.b);
	vector<double> b_ist(mesh.num_nod);
	solver.A.Ax(q, b_ist);
	for (int i = 0; i < mesh.num_nod; i++)
	{
		if (abs(slau.b[i] - b_ist[i]) > 1e-5)
			cout << " i " << i << " b[i] " << slau.b[i] << " Aq[i] " << b_ist[i] << endl;
	}

	vector<double> q_ist(mesh.num_nod);
	for (int i = 0; i < mesh.num_nod; i++)
	{
		q_ist[i] = u_test(mesh.all_nodes[i].x, mesh.all_nodes[i].y, mesh.all_nodes[i].z);
	}
	solver.A.Ax(q_ist, b_ist);
	for (int i = 0; i < mesh.num_nod; i++)
	{
		if (abs(slau.b[i] - b_ist[i]) > 1e-5)
			cout << " i " << i << " b[i] " << slau.b[i] << " Aq[i] " << b_ist[i] << endl;
	}*/
	
	double residual = 0;
	for (int i = 0; i < mesh.num_nod; i++)
	{
		residual += abs(q[i] - u_test(mesh.all_nodes[i].x, mesh.all_nodes[i].y, mesh.all_nodes[i].z));
	}
	residual /= mesh.num_nod;
	cout << " error: " << residual << endl;

	return 0;
}

int FEM::GetG_Loc_test(int el_id) // получение локальной G
{
	double lambda = mesh.all_material[mesh.all_elems[el_id].obl].lambda;
	double
		hx = mesh.all_nodes[mesh.all_elems[el_id].node_loc[1]].x - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].x,
		hy = mesh.all_nodes[mesh.all_elems[el_id].node_loc[2]].y - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].y,
		hz = mesh.all_nodes[mesh.all_elems[el_id].node_loc[4]].z - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].z;

	vector<Node3D> nodes_coord(8);
	for (int i = 0; i < 8; i++)
	{
		nodes_coord[i].x = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].x;
		nodes_coord[i].y = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].y;
		nodes_coord[i].z = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].z;
	}
	gauss3method_3D.GetLocG(nodes_coord, A_loc);
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
			A_loc[i][j] *= lambda;
	}

	/*vector<vector<double>> G_loc(8, vector<double>(8));
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			G_loc[i][j] = lambda * (hy * hz / hx * G_linear[mu(i)][mu(j)] * M_linear[nu(i)][nu(j)] * M_linear[teta(i)][teta(j)] +
				hx * hz / hy * M_linear[mu(i)][mu(j)] * G_linear[nu(i)][nu(j)] * M_linear[teta(i)][teta(j)] +
				hy * hx / hz * M_linear[mu(i)][mu(j)] * M_linear[nu(i)][nu(j)] * G_linear[teta(i)][teta(j)]);
		}
	}
	double diffrence = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			diffrence += abs(A_loc[i][j] - G_loc[i][j]);
		}
	}
	if (diffrence / 16 > 1e-14)
	{
		cout << " !!! warring. on elem " << el_id << " diffrence in G is " << diffrence / 16 << endl;
	}*/

	return 0;

}

int FEM::Getb_Loc_test(int el_id) // получение локального b
{
	for (int i = 0; i < 8; i++)
	{
		b_loc[i] = 0;
	}
	double x1, x2, y1, y2, z1, z2;
	x1 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].x;
	x2 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[1]].x;
	y1 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].y;
	y2 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[2]].y;
	z1 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].z;
	z2 = mesh.all_nodes[mesh.all_elems[el_id].node_loc[4]].z;
	double f[] = {
		f_test(x1, y1, z1),
		f_test(x2, y1, z1),
		f_test(x1, y2, z1),
		f_test(x2, y2, z1),
		f_test(x1, y1, z2),
		f_test(x2, y1, z2),
		f_test(x1, y2, z2),
		f_test(x2, y2, z2) }; // верно

	vector<vector<double>> M_loc(8, vector<double>(8));
	vector<Node3D> nodes_coord(8);
	for (int i = 0; i < 8; i++)
	{
		nodes_coord[i].x = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].x;
		nodes_coord[i].y = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].y;
		nodes_coord[i].z = mesh.all_nodes[mesh.all_elems[el_id].node_loc[i]].z;
	}
	gauss3method_3D.GetLocM(nodes_coord, M_loc);

	/*double
		hx = mesh.all_nodes[mesh.all_elems[el_id].node_loc[1]].x - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].x,
		hy = mesh.all_nodes[mesh.all_elems[el_id].node_loc[2]].y - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].y,
		hz = mesh.all_nodes[mesh.all_elems[el_id].node_loc[4]].z - mesh.all_nodes[mesh.all_elems[el_id].node_loc[0]].z;
	vector<vector<double>> test_M_loc(8, vector<double>(8));

	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			test_M_loc[i][j] = hx * hz * hy * M_linear[mu(i)][mu(j)] * M_linear[nu(i)][nu(j)] * M_linear[teta(i)][teta(j)];
		}
	}
	double diffrence = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			diffrence += abs(M_loc[i][j] - test_M_loc[i][j]);
		}
	}
	if (diffrence / 16 > 1e-14)
	{
		cout << " !!! warring. on elem " << el_id << " diffrence in M is " << diffrence / 16 << endl;
	}*/


	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			b_loc[i] += f[j] * M_loc[i][j];
		}
	}

	double gamma = mesh.all_material[mesh.all_elems[el_id].obl].lambda;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			A_loc[i][j] += gamma * M_loc[i][j];
		}
	}
	return 0;
}

int FEM::WriteP(string path)
{
	ofstream file;
	file.open(path + "result.txt");
	for (int i = 0; i < q.size(); i++)
	{
		file << q[i] << endl;
	}
	file.close();
	file.clear();

	return 0;
}