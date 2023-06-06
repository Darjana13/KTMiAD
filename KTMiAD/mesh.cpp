#include "mesh.h"

vector<set<int>> face2node =
{
	{0, 2, 4, 6},
	{1, 3, 5, 7},
	{0, 1, 4, 5},
	{2, 3, 6, 7},
	{0, 1, 2, 3},
	{4, 5, 6, 7},
};

vector<map<int, int>> node2edge = // node2edge[node_loc_max][node_loc_min] = edge_loc
{
	{},
	{{0, 0}},
	{{0, 1}},
	{{1, 2}, {2, 3}},
	{{0, 4}},
	{{4, 8},{1, 5}},
	{{4, 9},{2, 6}},
	{{3, 7},{5, 10},{6, 11}}
};
vector<vector<int>> edge2node =
{
	{0, 1},
	{0, 2},
	{1, 3},
	{2, 3},
	{0, 4},
	{1, 5},
	{2, 6},
	{3, 7},
	{4, 5},
	{4, 6},
	{5, 7},
	{6, 7}

};

vector<vector<int>> bf_id =
{
	{0, 2, 6, 8, 18, 20, 24, 26}, // node
	{12, 14, 10, 16, 4, 22}, // fc
	{1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25}, // eg
	{13}  // elem
};

int mesh::ReadMat()
{
	all_material.resize(1);
	all_material[0].lambda = 1;
	all_material[0].gamma = 1;

	return 0;
}
int mesh::ReadMesh()
{
	ReadMat();

	int kx, ky, kz, No;

	ifstream file;
	file.open("points.txt");

	file >> kx >> ky;

	x_base_line.resize(kx);
	for (int i = 0; i < kx; i++)
	{
		x_base_line[i].points.resize(ky);
	}

	y_base_line.resize(ky);
	for (int i = 0; i < ky; i++)
	{
		y_base_line[i].points.resize(kx);
	}

	for (int i = 0; i < ky; i++)
	{
		for (int j = 0; j < kx; j++)
		{
			file >> x_base_line[j].points[i] >> y_base_line[i].points[j];
		}
	}

	file >> kz;
	z_base_line.resize(kz);

	for (int i = 0; i < kz; i++)
	{
		z_base_line[i].points.resize(1);

		file >> z_base_line[i].points[0];
	}

	// считываем материалы областей
	file >> No;
	area.resize(No);
	for (int i = 0; i < No; i++)
	{
		file >> area[i].material >> area[i].x_start >> area[i].x_end >> area[i].y_start >> area[i].y_end >> area[i].z_start >> area[i].z_end;
	}

	// нужно считать кривые линии
	// считать число кривых
	// первое число это 1 или -1 показывает вертикальное или горизонтальне направление
	// потом 2 числа это номер ломаной и номер ее участка
	// последние 2 числа это центр окружности по х и у

	// разбиение сетки 
	// по идее тут три строки для x y z
	// возможно, нужно идти по основным узам, определяя участок между ними
	// и для каждого рассматриваемого участка считаем на сколько интервалов его разбить и с каким коэф
	vector<int> offset_x(kx), offset_y(ky), offset_z(kz);

	coef_x.resize(kx);
	offset_x[0] = 0;
	for (int i = 0; i < kx - 1; i++)
	{
		file >> coef_x[i].n >> coef_x[i].k;
		offset_x[i + 1] = offset_x[i] + coef_x[i].n;
		coef_x[i].start = offset_x[i];
		coef_x[i].end = offset_x[i + 1];
	}

	coef_y.resize(ky);
	offset_y[0] = 0;
	for (int i = 0; i < ky - 1; i++)
	{
		file >> coef_y[i].n >> coef_y[i].k;
		offset_y[i + 1] = offset_y[i] + coef_y[i].n;
		coef_y[i].start = offset_y[i];
		coef_y[i].end = offset_y[i + 1];
	}

	coef_z.resize(kz);
	offset_z[0] = 0;
	for (int i = 0; i < kz - 1; i++)
	{
		file >> coef_z[i].n >> coef_z[i].k;
		offset_z[i + 1] = offset_z[i] + coef_z[i].n;
		coef_z[i].start = offset_z[i];
		coef_z[i].end = offset_z[i + 1];
	}
	for (int i = 0; i < No; i++)
	{
		area[i].x_start = offset_x[area[i].x_start];
		area[i].x_end = offset_x[area[i].x_end];
		area[i].y_start = offset_y[area[i].y_start];
		area[i].y_end = offset_y[area[i].y_end];
		area[i].z_start = offset_z[area[i].z_start];
		area[i].z_end = offset_z[area[i].z_end];
	}

	x_line.resize(offset_x[kx - 1] + 1);
	for (int i = 0; i < offset_x[kx - 1] + 1; i++)
		x_line[i].points.resize(offset_y[ky - 1] + 1);
	y_line.resize(offset_y[ky - 1] + 1);
	for (int i = 0; i < offset_y[ky - 1] + 1; i++)
		y_line[i].points.resize(offset_x[kx - 1] + 1);
	z_line.resize(offset_z[kz - 1] + 1);
	for (int i = 0; i < offset_z[kz - 1] + 1; i++)
		z_line[i].points.resize(1);

	c_line tmp;
	type_coef_lines tmp_coef;
	// дробим y_line по x
	for (int i = 0; i < ky; i++) 
	{
		for (int j = 0; j < kx - 1; j++)
		{
			y_line[offset_y[i]].points[offset_x[j]] = y_base_line[i].points[j];
			y_line[offset_y[i]].points[offset_x[j + 1]] = y_base_line[i].points[j + 1];
			SplitSegment(coef_x[j], y_line[offset_y[i]]);
		}
		if (i != 0)
		{
			tmp.points.resize(coef_y[i - 1].n + 1);
			tmp_coef = coef_y[i - 1];
			tmp_coef.start = 0;
			tmp_coef.end = coef_y[i - 1].n;

			for (int j = 0; j < x_line.size(); j++)
			{
				tmp.points[0] = y_line[offset_y[i - 1]].points[j];
				tmp.points[coef_y[i - 1].n] = y_line[offset_y[i]].points[j];
				SplitSegment(tmp_coef, tmp);
				for (int l = 0; l < coef_y[i - 1].n; l++)
				{
					y_line[offset_y[i - 1] + l].points[j] = tmp.points[l];
				}
			}
		}
	}
	// дробим x_line по y
	for (int i = 0; i < kx; i++)
	{
		for (int k = 0; k < ky - 1; k++)
		{
			x_line[offset_x[i]].points[offset_y[k]] = x_base_line[i].points[k];
			x_line[offset_x[i]].points[offset_y[k+1]] = x_base_line[i].points[k + 1];
			SplitSegment(coef_y[k], x_line[offset_x[i]]);
		}
		if (i != 0)
		{
			tmp.points.resize(coef_x[i - 1].n + 1);
			tmp_coef = coef_x[i - 1];
			tmp_coef.start = 0;
			tmp_coef.end = coef_x[i - 1].n;

			for (int j = 0; j < y_line.size(); j++)
			{
				tmp.points[0] = x_line[offset_x[i - 1]].points[j];
				tmp.points[coef_x[i - 1].n] = x_line[offset_x[i]].points[j];
				SplitSegment(tmp_coef, tmp);
				for (int l = 0; l < coef_x[i - 1].n; l++)
				{
					x_line[offset_x[i - 1] + l].points[j] = tmp.points[l];
				}
			}
		}
	}



	for (int i = 0; i < kz - 1; i++)
	{
		if (coef_z[i].k == 1)
		{
			double step;
			step = z_base_line[i + 1].points[0] - z_base_line[i].points[0];
			step = step / coef_z[i].n;
			for (int j = 0; j < coef_z[i].n; j++)
			{
				z_line[offset_z[i] + j].points[0] = z_base_line[i].points[0] + step * j;
			}
		}
		else
		{
			double step;
			step = z_base_line[i + 1].points[0] - z_base_line[i].points[0];
			step = step * (coef_z[i].k - 1) / (pow(coef_z[i].k, coef_z[i].n) - 1);
			z_line[offset_z[i]].points[0] = z_base_line[i].points[0];
			for (int p = 0; p < coef_z[i].n; p++)
			{
				z_line[offset_z[i] + p + 1].points[0] = z_line[offset_z[i] + p].points[0] + step * pow(coef_z[i].k, p);
			}
		}
		z_line[offset_z[kz - 1]].points[0] = z_base_line[kz - 1].points[0];
	}
	//________________________ в этом месте считаем количество узлов и элементов _____________________________
	num_nod = x_line.size() * y_line.size() * z_line.size();
	num_el = (x_line.size() - 1) * (y_line.size() - 1) * (z_line.size() - 1);

	return 0;
}

// lines[s_id] && lines[e_id] exist
int mesh::SplitSegment(const type_coef_lines& coef, c_line& lines)
{
	if (coef.k == 1)
	{
		double step;
		step = lines.points[coef.end] - lines.points[coef.start];
		step = step / coef.n;
		for (int i = 1; i < coef.n; i++)
		{
			lines.points[coef.start + i] = lines.points[coef.start] + step * i;
		}
	}
	else
	{
		double step;
		step = lines.points[coef.end] - lines.points[coef.start];
		step = step * (coef.k - 1) / (pow(coef.k, coef.n) - 1);
		//for (int p = 0; p < coef.n - 1; p++)
		for (int p = 0; p < coef.n; p++)
		{
			lines.points[coef.start + p + 1] = lines.points[coef.start + p] + step * pow(coef.k, p);
		}
	}

	return 0;
}

int mesh::InitMemory() // пункт первый последняя часть выделение памяти под массив номеров элементов и массив локальных узлов 
{
	all_elems.resize(num_el);
	for (int i = 0; i < num_el; i++)
	{
		all_elems[i].node_loc.resize(8);
		all_elems[i].face_loc.resize(6);
		all_elems[i].edge_loc.resize(12);
		all_elems[i].bf_loc.resize(27);

		all_elems[i].obl = -1;

	}
	z_size = z_line.size();
	y_size = y_line.size();
	x_size = x_line.size();
	all_nodes.reserve(num_nod);
	face_to_elem.reserve(6 * num_el);
	edge_to_elem.reserve(12 * num_el);
	return 0;
}



int type_coef_lines::SplitSegment(vector<c_line>& base, vector<c_line>& result, int size_points)
{
	int cur_i = result.size();
	result.resize(cur_i + n);
	for (int i = 0; i < cur_i + n - 1; i++)
		result[i].points.resize(size_points);
	result[cur_i + n - 1] = base[end];

	if (k == 1)
	{
		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
		{
			double step;
			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
			step = step / n;
			for (int i = 1; i < n; i++)
			{
				result[cur_i - 1 + i].points[j] = result[cur_i - 1].points[j] + step * i;
			}
		}
	}
	else
	{
		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
		{
			double step;
			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
			step = step * (k - 1) / (pow(k, n) - 1);
			for (int p = 0; p < n - 1; p++)
			{
				result[cur_i + p].points[j] = result[cur_i + p - 1].points[j] + step * pow(k, p);

			}
		}
	}

	return 0;
}

//int type_coef_lines::SplitSegment(vector<c_line>& base, vector<c_line>& result)
//{
//	int cur_i = result.size();
//	result.resize(cur_i + n);
//	for (int i = 0; i < cur_i + n - 1; i++)
//		result[i].points.resize(2);
//	result[cur_i + n - 1] = base[end];
//
//	if (k == 1)
//	{
//		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
//		{
//			double step;
//			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
//			step = step / n;
//			for (int i = 1; i < n; i++)
//			{
//				result[cur_i - 1 + i].points[j] = result[cur_i - 1].points[j] + step * i;
//			}
//		}
//	}
//	else
//	{
//		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
//		{
//			double step;
//			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
//			step = step * (k - 1) / (pow(k, n) - 1);
//			for (int p = 0; p < n - 1; p++)
//			{
//				result[cur_i + p].points[j] = result[cur_i + p - 1].points[j] + step * pow(k, p);
//
//			}
//		}
//	}
//
//	return 0;
//}
//int type_coef_lines::SplitSegment_z(vector<c_line>& base, vector<c_line>& result)
//{
//	int cur_i = result.size();
//	result.resize(cur_i + n);
//	for (int i = 0; i < cur_i + n - 1; i++)
//		result[i].points.resize(1);
//	result[cur_i + n - 1] = base[end];
//
//	if (k == 1)
//	{
//		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
//		{
//			double step;
//			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
//			step = step / n;
//			for (int i = 1; i < n; i++)
//			{
//				result[cur_i - 1 + i].points[j] = result[cur_i - 1].points[j] + step * i;
//			}
//		}
//	}
//	else
//	{
//		for (int j = 0; j < result[cur_i - 1].points.size(); j++)
//		{
//			double step;
//			step = result[cur_i + n - 1].points[j] - result[cur_i - 1].points[j];
//			step = step * (k - 1) / (pow(k, n) - 1);
//			for (int p = 0; p < n - 1; p++)
//			{
//				result[cur_i + p].points[j] = result[cur_i + p - 1].points[j] + step * pow(k, p);
//
//			}
//		}
//	}
//
//	return 0;
//}


int mesh::FillArea()
{
	int z_size = z_line.size();
	int y_size = y_line.size();
	int x_size = x_line.size();

	int k_area = area.size();
	int i_elem = 0;

	for (int i_area = 0; i_area < k_area; i_area++)
	{
		for (int z = area[i_area].z_start; z < area[i_area].z_end; z++)
			for (int y = area[i_area].y_start; y < area[i_area].y_end; y++)
				for (int x = area[i_area].x_start; x < area[i_area].x_end; x++)
				{
					int i = 0;
					i = x + y * x_size + z * x_size * y_size;
					all_elems[i_elem].obl = i_area;
					all_elems[i_elem].node_loc[0] = i;
					all_elems[i_elem].node_loc[1] = i + 1;
					all_elems[i_elem].node_loc[2] = i + x_size;
					all_elems[i_elem].node_loc[3] = i + x_size + 1;
					all_elems[i_elem].node_loc[4] = i + x_size * y_size;
					all_elems[i_elem].node_loc[5] = i + x_size * y_size + 1;
					all_elems[i_elem].node_loc[6] = i + x_size * y_size + x_size;
					all_elems[i_elem].node_loc[7] = i + x_size * y_size + x_size + 1;
					i_elem++;
				}
	}



	return 0;
}

int mesh::AddNodeToRealNodes(int glob_num, int i_node)
{
	type_points p;
	int real_num;
	if (glob2real.count(glob_num))
		real_num = glob2real[glob_num];
	else
	{
		glob2real[glob_num] = i_node;
		p.x = x_line[(glob_num % (x_size * y_size)) % x_size].points[(glob_num % (x_size * y_size)) / x_size];
		p.y = y_line[(glob_num % (x_size * y_size)) / x_size].points[(glob_num % (x_size * y_size)) % x_size];
		p.z = z_line[glob_num / (x_size * y_size)].points[0];
		all_nodes.push_back(p);
		real_num = i_node;
		//i_node++;
	}
	return real_num;
}

int mesh::FillArea2()
{
	int x_eg_all = (x_size - 1) * (y_size) * (z_size);
	int	y_eg_all = (y_size - 1) * (x_size) * (z_size);
	int x_fc_all = x_size * (y_size - 1) * (z_size - 1);
	int	y_fc_all = y_size * (x_size - 1) * (z_size - 1);

	int k_area = area.size();
	int i_elem = 0;
	int i_node = 0;
	int glob_num_node, glob_num_face, glob_num_edge, glob_num_bf;
	int real_num_node, real_num_face, real_num_edge, real_num_bf;
	type_points p;

	int i = 0, i_fc = 0, i_edge = 0, i_bf = 0;

	for (int i_area = 0; i_area < k_area; i_area++)
	{
		for (int z = area[i_area].z_start; z < area[i_area].z_end; z++)
			for (int y = area[i_area].y_start; y < area[i_area].y_end; y++)
				for (int x = area[i_area].x_start; x < area[i_area].x_end; x++)
				{
					i = x + y * x_size + z * x_size * y_size;
					all_elems[i_elem].obl = area[i_area].material;
#pragma region NumNodes
					//all_elems[i_elem].node_loc[0] = i;
					glob_num_node = i;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[0] = real_num_node;

					//all_elems[i_elem].node_loc[1] = i + 1;
					glob_num_node = i + 1;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[1] = real_num_node;

					// all_elems[i_elem].node_loc[2] = i + x_size;
					glob_num_node = i + x_size;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[2] = real_num_node;

					//   all_elems[i_elem].node_loc[3] = i + x_size + 1;
					glob_num_node = i + x_size + 1;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[3] = real_num_node;

					// all_elems[i_elem].node_loc[4] = i + x_size * y_size;
					glob_num_node = i + x_size * y_size;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[4] = real_num_node;

					//all_elems[i_elem].node_loc[5] = i + x_size * y_size + 1;
					glob_num_node = i + x_size * y_size + 1;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[5] = real_num_node;

					//all_elems[i_elem].node_loc[6] = i + x_size * y_size + x_size;
					glob_num_node = i + x_size * y_size + x_size;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[6] = real_num_node;

					//all_elems[i_elem].node_loc[7] = i + x_size * y_size + x_size + 1;
					glob_num_node = i + x_size * y_size + x_size + 1;
					real_num_node = AddNodeToRealNodes(glob_num_node, i_node);
					if (real_num_node == i_node) i_node++;
					all_elems[i_elem].node_loc[7] = real_num_node;
#pragma endregion
#pragma region NumFaces
					// x_fc_all = x_size * (y_size - 1) * (z_size - 1)
					// y_fc_all = y_size * (x_size - 1) * (z_size - 1)
					//all_elems[i_elem].face_loc[0] = x + y * (x_size) + z * (x_size) * (y_size - 1);
					//all_elems[i_elem].face_loc[1] = all_elems[i_elem].face_loc[0] + 1;
					//all_elems[i_elem].face_loc[2] = x_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size);
					//all_elems[i_elem].face_loc[3] = all_elems[i_elem].face_loc[2] + (x_size - 1);
					//all_elems[i_elem].face_loc[4] = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1);
					//all_elems[i_elem].face_loc[5] = all_elems[i_elem].face_loc[4] + (x_size - 1);

					//all_elems[i_elem].face_loc[0] =  x + y * (x_size) + z * (x_size) * (y_size - 1);
					glob_num_face = x + y * (x_size)+z * (x_size) * (y_size - 1);
					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[0] = real_num_face;
					//all_elems[i_elem].face_loc[1] = all_elems[i_elem].face_loc[0] + 1;
					glob_num_face = x + y * (x_size)+z * (x_size) * (y_size - 1) + 1;
					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[1] = real_num_face;
					//all_elems[i_elem].face_loc[2] =x_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size);
					glob_num_face = x_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size);
					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[2] = real_num_face;
					//all_elems[i_elem].face_loc[3] = all_elems[i_elem].face_loc[2] + (x_size - 1);
					glob_num_face = x_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size)+(x_size - 1);
					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[3] = real_num_face;
					//all_elems[i_elem].face_loc[4] = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1);
					//glob_num_face = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1);
					glob_num_face = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1);
					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[4] = real_num_face;
					//all_elems[i_elem].face_loc[5] = all_elems[i_elem].face_loc[4] + (x_size - 1);
					//glob_num_face = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1) + (x_size - 1);
					glob_num_face = x_fc_all + y_fc_all + x + y * (x_size - 1) + z * (x_size - 1) * (y_size - 1) + (x_size - 1) * (y_size - 1);

					if (face_glob2real.count(glob_num_face))
					{
						real_num_face = face_glob2real[glob_num_face];
						face_to_elem[real_num_face].insert(i_elem);
					}
					else
					{
						face_glob2real[glob_num_face] = i_fc;
						real_num_face = i_fc;
						i_fc++;
						set<int> neib = { i_elem };
						face_to_elem.push_back(neib);
					}
					all_elems[i_elem].face_loc[5] = real_num_face;

#pragma endregion
#pragma region NumEdges
					// x_eg_all = (x_size - 1) * (y_size) * (z_size)
					// y_eg_all = (y_size - 1) * (x_size) * (z_size)
					//all_elems[i_elem].edge_loc[0] = x + y * (x_size-1) + z * (x_size-1) * (y_size);
					//all_elems[i_elem].edge_loc[1] = x_eg_all + x + y * (x_size) + z * (y_size-1) * (x_size) ;
					//all_elems[i_elem].edge_loc[2] = all_elems[i_elem].edge_loc[1] + 1;
					//all_elems[i_elem].edge_loc[3] = all_elems[i_elem].edge_loc[0] + (x_size-1);
					//all_elems[i_elem].edge_loc[4] = x_eg_all + y_eg_all + x + y *x_size + z * x_size * y_size ;
					//all_elems[i_elem].edge_loc[5] = all_elems[i_elem].edge_loc[4]  + 1;
					//all_elems[i_elem].edge_loc[6] = all_elems[i_elem].edge_loc[4] + x_size ;
					//all_elems[i_elem].edge_loc[7] = all_elems[i_elem].edge_loc[4] + x_size + 1;
					//all_elems[i_elem].edge_loc[8] = all_elems[i_elem].edge_loc[0] + (x_size-1) * (y_size);
					//all_elems[i_elem].edge_loc[9] = all_elems[i_elem].edge_loc[1] + (y_size-1) * (x_size);
					//all_elems[i_elem].edge_loc[10] = all_elems[i_elem].edge_loc[9] + 1;
					//all_elems[i_elem].edge_loc[11] = all_elems[i_elem].edge_loc[8] +  (x_size-1);

					//all_elems[i_elem].edge_loc[0] = x + y * (x_size-1) + z * (x_size-1) * (y_size);
					glob_num_edge = x + y * (x_size - 1) + z * (x_size - 1) * (y_size);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[0] = real_num_edge;

					//all_elems[i_elem].edge_loc[1] = x_eg_all + y + x * (y_size-1) + z * (y_size-1) * (x_size) ;
					glob_num_edge = x_eg_all + y + x * (y_size - 1) + z * (y_size - 1) * (x_size);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[1] = real_num_edge;
					//all_elems[i_elem].edge_loc[2] = all_elems[i_elem].edge_loc[1] + (y_size-1);
					glob_num_edge = x_eg_all + y + x * (y_size - 1) + z * (y_size - 1) * (x_size)+(y_size - 1);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[2] = real_num_edge;
					//all_elems[i_elem].edge_loc[3] = all_elems[i_elem].edge_loc[0] + (x_size-1);
					glob_num_edge = x + y * (x_size - 1) + z * (x_size - 1) * (y_size)+(x_size - 1);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[3] = real_num_edge;
					
					//all_elems[i_elem].edge_loc[4] = x_eg_all + y_eg_all + x + y *x_size + z * x_size * y_size ;
					glob_num_edge = x_eg_all + y_eg_all + x + y * x_size + z * x_size * y_size;
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[4] = real_num_edge;
					//all_elems[i_elem].edge_loc[5] = all_elems[i_elem].edge_loc[4]  + 1;
					glob_num_edge = x_eg_all + y_eg_all + x + y * x_size + z * x_size * y_size + 1;
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[5] = real_num_edge;
					//all_elems[i_elem].edge_loc[6] = all_elems[i_elem].edge_loc[4] + x_size ;
					glob_num_edge = x_eg_all + y_eg_all + x + y * x_size + z * x_size * y_size + x_size;
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[6] = real_num_edge;
					//all_elems[i_elem].edge_loc[7] = all_elems[i_elem].edge_loc[4] + x_size + 1;
					glob_num_edge = x_eg_all + y_eg_all + x + y * x_size + z * x_size * y_size + x_size + 1;
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[7] = real_num_edge;
					//all_elems[i_elem].edge_loc[8] = all_elems[i_elem].edge_loc[0] + (x_size-1) * (y_size);
					glob_num_edge = x + y * (x_size - 1) + z * (x_size - 1) * (y_size)+(x_size - 1) * (y_size);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[8] = real_num_edge;
					//all_elems[i_elem].edge_loc[9] = all_elems[i_elem].edge_loc[1] + (y_size-1) * (x_size);
					glob_num_edge = x_eg_all + y + x * (y_size - 1) + z * (y_size - 1) * (x_size)+(y_size - 1) * (x_size);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[9] = real_num_edge;
					//all_elems[i_elem].edge_loc[10] = all_elems[i_elem].edge_loc[9] + (y_size-1);
					glob_num_edge = x_eg_all + y + x * (y_size - 1) + z * (y_size - 1) * (x_size)+(y_size - 1) * (x_size)+(y_size - 1);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[10] = real_num_edge;
					//all_elems[i_elem].edge_loc[11] = all_elems[i_elem].edge_loc[8] +  (x_size-1);
					glob_num_edge = x + y * (x_size - 1) + z * (x_size - 1) * (y_size)+(x_size - 1) * (y_size)+(x_size - 1);
					if (edge_glob2real.count(glob_num_edge))
					{
						real_num_edge = edge_glob2real[glob_num_edge];
						edge_to_elem[real_num_edge].insert(i_elem);
					}
					else
					{
						edge_glob2real[glob_num_edge] = i_edge;
						real_num_edge = i_edge;
						set<int> tmp = { i_elem };
						edge_to_elem.push_back(tmp);
						i_edge++;
					}
					all_elems[i_elem].edge_loc[11] = real_num_edge;

#pragma endregion
					i_elem++;
				}
	}

	num_el = i_elem;
	num_nod = i_node;
	eg_all = i_edge;
	fc_all = i_fc;

	
	i_bf = 0;
	for (int el_id = 0; el_id < num_el; el_id++)
	{
#pragma region NumBF
		// node = 0...num_nod - 1
		// fc = num_nod...fc_all - 1
		// edge = num_node + fc_all...eg_all - 1
		// el = num_node + fc_all + eg_all
		for (i = 0; i < 8; i++)
		{
			glob_num_bf = all_elems[el_id].node_loc[i];
			if (bf_glob2real.count(glob_num_bf))
			{
				real_num_bf = bf_glob2real[glob_num_bf];
			}
			else
			{
				bf_glob2real[glob_num_bf] = i_bf;
				real_num_bf = i_bf;
				i_bf++;
			}
			all_elems[el_id].bf_loc[bf_id[0][i]] = real_num_bf;
		}
		for (i = 0; i < 6; i++)
		{
			glob_num_bf = num_nod + all_elems[el_id].face_loc[i];
			if (bf_glob2real.count(glob_num_bf))
			{
				real_num_bf = bf_glob2real[glob_num_bf];
			}
			else
			{
				bf_glob2real[glob_num_bf] = i_bf;
				real_num_bf = i_bf;
				i_bf++;
			}
			all_elems[el_id].bf_loc[bf_id[1][i]] = real_num_bf;
		}
		for (i = 0; i < 12; i++)
		{
			glob_num_bf = num_nod + fc_all + all_elems[el_id].edge_loc[i];
			if (bf_glob2real.count(glob_num_bf))
			{
				real_num_bf = bf_glob2real[glob_num_bf];
			}
			else
			{
				bf_glob2real[glob_num_bf] = i_bf;
				real_num_bf = i_bf;
				i_bf++;
			}
			all_elems[el_id].bf_loc[bf_id[2][i]] = real_num_bf;
		}
		glob_num_bf = num_nod + fc_all + eg_all + el_id;
		if (bf_glob2real.count(glob_num_bf))
		{
			real_num_bf = bf_glob2real[glob_num_bf];
		}
		else
		{
			bf_glob2real[glob_num_bf] = i_bf;
			real_num_bf = i_bf;
			i_bf++;
		}
		all_elems[el_id].bf_loc[bf_id[3][0]] = real_num_bf;
#pragma endregion
	}
	bf_all = i_bf;

	return 0;
}

int mesh::GetEdge(int p1, int p2)
{
	int flag1 = -1, flag2 = -1;
	int find_elem;
	for (find_elem = 0; find_elem < num_el && (flag1 == -1 && flag2 == -1); find_elem++)
	{
		flag1 = -1;
		flag2 = -1;
		for (int j = 0; j < 8 && (flag1 == -1 || flag2 == -1); j++)
		{
			if (all_elems[find_elem].node_loc[j] == p1)
				flag1 = j;
			else if (all_elems[find_elem].node_loc[j] == p2)
				flag2 = j;
		}
	}
	if (flag1 == -1 && flag2 == -1)
	{
		cout << "error: no edge (can't find elem)" << endl;
		return -1;
	}
	find_elem--;
	
	if (flag1 < flag2)
		return all_elems[find_elem].edge_loc[node2edge[flag2][flag1]];
	else
		return all_elems[find_elem].edge_loc[node2edge[flag1][flag2]];
}

int mesh::GetNeibEdge(int edge, int &p1, int& p2, set<int> &elems)
{
	elems = edge_to_elem[edge];
	int cur_el = *(edge_to_elem[edge].begin());
	for (int i = 0; i < 12; i++)
	{
		if (all_elems[cur_el].edge_loc[i] == edge)
		{
			p1 = edge2node[i][0];
			p2 = edge2node[i][1];
			return 0;
		}
	}
	cout << "error: no edge " << endl;
	return 1;
}

int mesh::GetFace(int p1, int p2, int p3, int p4)
{
	int flag1 = -1, flag2 = -1, flag3 = -1, flag4 = -1;
	int find_elem;
	for (find_elem = 0; find_elem < num_el && (flag1 == -1 && flag2 == -1 && flag3 == -1 && flag4 == -1); find_elem++)
	{
		flag1 = -1;
		flag2 = -1;
		flag3 = -1;
		flag4 = -1;
		for (int j = 0; j < 8 && (flag1 == -1 || flag2 == -1 || flag3 == -1 || flag4 == -1); j++)
		{
			if (all_elems[find_elem].node_loc[j] == p1)
				flag1 = j;
			else if (all_elems[find_elem].node_loc[j] == p2)
				flag2 = j;
			else if (all_elems[find_elem].node_loc[j] == p3)
				flag3 = j;
			else if (all_elems[find_elem].node_loc[j] == p4)
				flag4 = j;
		}
	}
	if (flag1 == -1 || flag2 == -1 || flag3 == -1 || flag4 == -1)
	{
		cout << "error: no face (can't find elem)" << endl;
		return -1;
	}
	find_elem--;

	set<int> node;
	node.insert(flag1);
	node.insert(flag2);
	node.insert(flag3);
	node.insert(flag4);
	for (int i = 0; i < 6; i++)
	{
		if (face2node[i] == node)
			return all_elems[find_elem].face_loc[i];
	}
	cout << "error: no face " << endl;
	return -1;
}

int mesh::GetNeibFace(int face, set<int>& nodes, set<int> &elems)
{
	elems = face_to_elem[face];
	int cur_el = *(face_to_elem[face].begin());
	for (int i = 0; i < 6; i++)
	{
		if (all_elems[cur_el].face_loc[i] == face)
		{
			auto iterator = face2node[i].begin();
			nodes.clear();

			for (; iterator != face2node[i].end(); iterator++)
			{
				nodes.insert(all_elems[cur_el].node_loc[(*iterator)]);
			}
			return 0;
		}
	}
	cout << "error: no face in neib " << endl;
	return 1;
}

//vector<vector<int>> bf_id =
//{
//	{0, 2, 6, 8, 18, 20, 24, 26}, // node
//	{12, 14, 10, 16, 4, 22}, // fc
//	{1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25}, // eg
//	{13}  // elem
//};
type_points mesh::GetCoordBF(int el, int bf_loc)
{
	type_points p;

	set<int> node = { 0, 2, 6, 8, 18, 20, 24, 26 };
	set<int> fc = { 12, 14, 10, 16, 4, 22 };
	set<int> eg = { 1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25 };
	// elem = {13}
	vector<int> bf2obj = { 0, 0, 1, 1, 4, 2, 2, 3, 3, 4, 2, 5, 0, 0, 1, 6, 3, 7, 4, 8, 5, 9, 5, 10, 6, 11, 7 };

	if (node.count(bf_loc))
	{
		return all_nodes[all_elems[el].node_loc[bf2obj[bf_loc]]];
	}
	else if (fc.count(bf_loc))
	{
		int fc_loc = bf2obj[bf_loc];
		p = GetCentreFace(el, fc_loc);
		return p;
	}
	else if (eg.count(bf_loc))
	{
		int eg_loc = bf2obj[bf_loc];
		int p1 = edge2node[eg_loc][0];
		int p2 = edge2node[eg_loc][1];
		p = all_nodes[all_elems[el].node_loc[p1]].GetPoint(all_nodes[all_elems[el].node_loc[p2]]);
		return p;
	}
	else
	{
		p = GetCentreElem(el);
		return p;
	}
}

type_points mesh::GetCentreFace(int el, int fc_loc)
{
	type_points p;
	p.x = 0;
	p.y = 0;
	p.z = 0;

	auto f2n = face2node[fc_loc].begin();
	for (; f2n != face2node[fc_loc].end(); f2n++)
	{
		p.x += all_nodes[all_elems[el].node_loc[(*f2n)]].x;
		p.y += all_nodes[all_elems[el].node_loc[(*f2n)]].y;
		p.z += all_nodes[all_elems[el].node_loc[(*f2n)]].z;

	}

	p.x /= 4;
	p.y /= 4;
	p.z /= 4;

	return p;
}

type_points mesh::GetCentreElem(int el)
{
	type_points p;
	p.x = 0;
	p.y = 0;
	p.z = 0;

	for (int i = 0; i < 8; i++)
	{
		p.x += all_nodes[all_elems[el].node_loc[i]].x;
		p.y += all_nodes[all_elems[el].node_loc[i]].y;
		p.z += all_nodes[all_elems[el].node_loc[i]].z;

	}

	p.x /= 8;
	p.y /= 8;
	p.z /= 8;

	return p;
}

int mesh::BuildS1()
{
	set<int> nodes, elems;
	for (int i = 0; i < fc_all; i++)
	{
		if (face_to_elem[i].size() == 1)
		{
			nodes.clear();
			GetNeibFace(i, nodes, elems);
			auto iterator = nodes.begin();
			for (; iterator != nodes.end(); iterator++)
			{
				S1.insert(make_pair((*iterator), 0));
			}
		}
	}


	return 0;
}

int mesh::WriteMesh(string path)
{
	ofstream file;

	file.open(path + "info.txt");
	file << "k_node " << num_nod << endl;
	file << "k_elem " << num_el << endl;
	file << "k_face " << fc_all << endl;
	file << "k_edge " << eg_all << endl;
	file << "k_bf   " << bf_all << endl;
	file << "k_S1   " << S1.size() << endl;

	file.close();
	file.clear();

	double
		x_min = 100000000000,
		x_max = -100000000000,
		y_min = 100000000000,
		y_max = -100000000000,
		z_min = 100000000000,
		z_max = -100000000000;
	file.open(path + "xyz_simple.txt");
	for (int i = 0; i < num_nod; i++)
	{
		file << all_nodes[i].x << " " << all_nodes[i].y << " " << all_nodes[i].z << endl;
		if (all_nodes[i].x < x_min)
			x_min = all_nodes[i].x;
		if (all_nodes[i].y < y_min)
			y_min = all_nodes[i].y;
		if (all_nodes[i].z < z_min)
			z_min = all_nodes[i].z;
		if (all_nodes[i].x > x_max)
			x_max = all_nodes[i].x;
		if (all_nodes[i].y > y_max)
			y_max = all_nodes[i].y;
		if (all_nodes[i].z > z_max)
			z_max = all_nodes[i].z;
	}
	file.close();
	file.clear();

	file.open(path + "xyz0_simple.txt");
	for (int i = 0; i < num_nod; i++)
	{
		// от -1 до 1
		file
			<< -1.0 + 2.0 * (all_nodes[i].x - x_min) / (x_max - x_min) << " "
			<< -1.0 + 2.0 * (all_nodes[i].y - y_min) / (y_max - y_min) << " "
			<< -1.0 + 2.0 * (all_nodes[i].z - z_min) / (z_max - z_min) << endl;
	}
	file.close();
	file.clear();

	file.open(path + "elem_simple.txt");
	for (int j = 0; j < num_el; j++)
	{
		for (int i = 0; i < 8; i++)
		{
			file << all_elems[j].node_loc[i] << " ";
		}
		file << all_elems[j].obl << endl;

	}
	file.close();
	file.clear();

	file.open(path + "elem_bf.txt");
	for (int j = 0; j < num_el; j++)
	{
		for (int i = 0; i < 27; i++)
		{
			file << all_elems[j].bf_loc[i] << " ";
		}
		file << endl;

	}
	file.close();
	file.clear();

	file.open(path + "elem_bf_with coord.txt");
	type_points p;
	for (int j = 0; j < num_el; j++)
	{
		file << " elem " << j << endl;
		for (int i = 0; i < 27; i++)
		{
			p = GetCoordBF(j, i);
			file << "      " << all_elems[j].bf_loc[i] << " x: " << p.x << " y: " << p.y << " z: " << p.z << endl;
		}
		file << endl;
	}
	file.close();
	file.clear();

	file.open(path + "elem_face.txt");
	for (int j = 0; j < num_el; j++)
	{
		for (int i = 0; i < 6; i++)
		{
			file << all_elems[j].face_loc[i] << " ";
		}
		file << endl;

	}
	file.close();
	file.clear();

	file.open(path + "elem_edge.txt");
	for (int j = 0; j < num_el; j++)
	{
		for (int i = 0; i < 12; i++)
		{
			file << all_elems[j].edge_loc[i] << " ";
		}
		file << endl;

	}
	file.close();
	file.clear();

	file.open(path + "face_neib.txt");
	for (int j = 0; j < fc_all; j++)
	{
		auto iterator = face_to_elem[j].begin();
		for (; iterator != face_to_elem[j].end(); iterator++)
		{
			file << (*iterator) << " ";
		}
		file << endl;
	}
	file.close();
	file.clear();

	return 0;
}