#include "mesh.h"

int SLAU::AddLocal(vector<vector<double>>& A_loc, vector<double>& b_loc, int el_id, mesh& mesh)
// внесение локальных A, b  в глобальную —Ћј”
{
	vector<int> L = mesh.all_elems[el_id].node_loc;
	int n_loc = mesh.all_elems[el_id].node_loc.size(); // размерность локальной матрицы
	int count = 0;

	for (int i = 0; i < n_loc; i++)
	{
		di[L[i]] += A_loc[i][i];
		b[L[i]] += b_loc[i];
		count += 1;
	}
	
	for (int i = 0; i < n_loc; i++)
	{
		int temp = ia[L[i]];
		//for (int j = 0; j < i; j++)
		for (int j = 0; j < n_loc; j++)
		{
			for (int k = temp; k < ia[L[i] + 1]; k++)
			{
				if (ja[k] == L[j])
				{
					al[k] += A_loc[j][i];
					au[k] += A_loc[i][j];
					count += 2;
					break;
				}
			}
		}
	}

	return 0;
}


int SLAU::SetS1_no_simmetry(mesh& mesh) // учет первых краевых
{
		for (auto iterator_S1 = mesh.S1.begin(); iterator_S1 != mesh.S1.end(); iterator_S1++)
		{
			int node_id = (*iterator_S1).first;
			di[node_id] = 1;
			b[node_id] = u_test(mesh.all_nodes[iterator_S1->first].x, mesh.all_nodes[iterator_S1->first].y, mesh.all_nodes[iterator_S1->first].z);
			for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
			{
				al[k] = 0; // строка, где 1ое краевое
			}
			for (int k = 0; k < ja.size(); k++)
			{
				if (ja[k] == node_id)
				{
					au[k] = 0; // строка, где 1ое краевое
				}
			}
		}
	return 0;
}

int SLAU::SetS1_test(mesh& mesh) // учет первых краевых
{
	int NS1 = mesh.S1.size();
	double max = 0;
	for (int i = 0; i < di.size(); i++)
	{
		if (di[i] > max)
			max = di[i];
	}
	for (auto iterator_S1 = mesh.S1.begin(); iterator_S1 != mesh.S1.end(); iterator_S1++)
	{
		di[iterator_S1->first] = max * BIG_VALUE;
		b[iterator_S1->first] = max * BIG_VALUE * u_test(mesh.all_nodes[iterator_S1->first].x, mesh.all_nodes[iterator_S1->first].y, mesh.all_nodes[iterator_S1->first].z);
	}
	return 0;
}

//int SLAU::GeneratePortret(mesh& mesh) // генераци€ портрета
//{
//	int N = mesh.num_nod;
//	int Kel = mesh.num_el;
//	int N_loc = 8;
//	ia.resize(N + 1);
//	ja.resize(N_loc* N_loc * Kel);
//	std::vector<int> temp_list1(N_loc* N_loc * Kel),
//		temp_list2(N_loc* N_loc * Kel);
//	std::vector<int> listbeg(N);
//	int listsize = 0;
//	for (int i = 0; i < N; i++)
//	{
//		listbeg[i] = 0;
//	}
//	for (int ielem = 0; ielem < Kel; ielem++)
//	{
//		for (int i = 0; i < N_loc; i++) // NumberOfUnknowns(ielem)?
//		{
//			int k = mesh.all_elems[ielem].node_loc[i];
//			for (int j = i + 1; j < N_loc; j++)// NumberOfUnknowns(ielem)?
//			{
//				int ind1 = k;
//				int ind2 = mesh.all_elems[ielem].node_loc[j];
//				if (ind2 < ind1)
//				{
//					ind1 = ind2;
//					ind2 = k;
//				}
//				int iaddr = listbeg[ind2];
//				if (iaddr == 0)
//				{
//					listsize++;
//					listbeg[ind2] = listsize;
//					temp_list1[listsize] = ind1;
//					temp_list2[listsize] = 0;
//				}
//				else
//				{
//					while (temp_list1[iaddr] < ind1 && temp_list2[iaddr] > 0)
//					{
//						iaddr = temp_list2[iaddr];
//					}
//					if (temp_list1[iaddr] > ind1)
//					{
//						listsize++;
//						temp_list1[listsize] = temp_list1[iaddr];
//						temp_list2[listsize] = temp_list2[iaddr];
//						temp_list1[iaddr] = ind1;
//						temp_list2[iaddr] = listsize;
//					}
//					else if (temp_list1[iaddr] < ind1)
//					{
//						listsize++;
//						temp_list2[iaddr] = listsize;
//						temp_list1[listsize] = ind1;
//						temp_list2[listsize] = 0;
//					}
//				}
//			}
//		}
//	}
//
//	ia[0] = 0;
//	for (int i = 0; i < N; i++)
//	{
//		ia[i + 1] = ia[i];
//		int iaddr = listbeg[i];
//		while (iaddr != 0)
//		{
//			ja[ia[i + 1]] = temp_list1[iaddr];
//			ia[i + 1]++;
//			iaddr = temp_list2[iaddr];
//		}
//	}
//
//	ja.resize(ia[N]);
//	au.resize(ja.size());
//	al.resize(ja.size());
//
//	di.resize(N);
//	b.resize(N);
//	return 0;
//}


int SLAU::GeneratePortret(mesh& mesh)
{
	int N = mesh.num_nod, Kel = mesh.num_el;
	int N_loc = mesh.all_elems[0].node_loc.size();
	ia.resize(N + 1);

	vector<set<int>> map(N);
	set<int> test;
	int k = 0;
	for (auto elem : mesh.all_elems)
	{
		for (auto i : elem.node_loc)
			for (auto j : elem.node_loc)
			{
				if (i == 7)
					test.insert(j);
				if (i > j)
					map[i].insert(j);
			}
	}
	ia[0] = 0;
	for (int i = 0; i < N; i++)
	{
		ia[i + 1] = ia[i] + map[i].size();
	}

	ja.resize(ia[N]);
	for (int i = 0; i < N; i++)
	{
		set <int> ::iterator it = map[i].begin();
		for (int j = 0; it != map[i].end(); it++, j++)
		{
			ja[ia[i] + j] = *(it);
		}
	}

	//a.resize(ja.size());
	au.resize(ja.size());
	al.resize(ja.size());

	di.resize(N);
	b.resize(N);

	return 0;
}

void SLAU::ClearValues()
{
	//fill(a.begin(), a.end(), 0.0);
	fill(au.begin(), au.end(), 0.0);
	fill(al.begin(), al.end(), 0.0);

	fill(b.begin(), b.end(), 0.0);
	fill(di.begin(), di.end(), 0.0);
}