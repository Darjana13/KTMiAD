#include "mesh.h"

__forceinline ostream& operator < (ostream& file, const double& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}
__forceinline ostream& operator < (ostream& file, const int& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file, double& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file, int& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

int main()
{
	FEM task;
	task.Init("");
	//task.mesh.WriteMesh("output\\");
	//return 0;

	task.SolveTask_test("");
	task.WriteP("output\\");
	task.mesh.WriteMesh("output\\");

   /*mesh test;
   test.ReadMesh();
   test.InitMemory();
   test.FillArea2();*/

   

   
   /*
   // взять номер ребра по номерам вершин
   int eg = test.GetEdge(1, 5);
   cout << " edge " << eg << endl;

   // по номеру ребра взять вершины и связанные с этим ребром элементы
   set<int> elems;
   int node_1, node_2;
   test.GetNeibEdge(eg, node_1, node_2, elems);
   auto iterator1 = elems.begin();
   cout << " nodes: " << node_1 << " " << node_2 << endl;
   cout << " elems: ";
   for (; iterator1 != elems.end(); iterator1++)
   {
	   cout << (*iterator1) << " ";
   }
   cout << endl;

   // взять номер грани по номерам вершин
   int fc = test.GetFace(1, 3, 5, 7);
   cout << " face " << fc << endl;

   // по номеру грани взять вершины и связанные с этой гранью элементы
   elems.clear();
   set<int> nodes;
   test.GetNeibFace(fc, nodes, elems);
   cout << " nodes: " ;
   iterator1 = nodes.begin();
   for (; iterator1 != nodes.end(); iterator1++)
   {
	   cout << (*iterator1) << " ";
   }
   cout << endl;

   iterator1 = elems.begin();
   cout << " elems: ";
   for (; iterator1 != elems.end(); iterator1++)
   {
	   cout << (*iterator1) << " ";
   }
   cout << endl;
   */

   return 0;
}