//
// Created by su on 2022/3/13.
//
#include<string>
#include<vector>
#include <map>

#include "basic_types.h"
#include "cross_field.h"
#include "basic_utils.h"

int readVTK(const std::string &name,
            std::vector<Vec3> &points,
            std::vector<ID2> &line,
            std::vector<ID3> &triangles,
            bool bigEndian = false) {
  std::vector<std::vector<std::vector<uint32_t>>> elements;

  //8种类别的元素，point，line，triangle，quad，tet，hex，prism，pyramid
  elements.resize(8);

  FILE *fp = fopen(name.c_str(), "rb");
  if (!fp) {
    fprintf(stdout, "Unable to open file '%s'\n", name.c_str());
    return 0;
  }

  char buffer[256], buffer2[256];

  if (!fgets(buffer, sizeof(buffer), fp)) {
    fclose(fp);
    return 0;
  } // version line
  if (!fgets(buffer, sizeof(buffer), fp)) {
    fclose(fp);
    return 0;
  } // title

  if (fscanf(fp, "%s", buffer) != 1) // ASCII or BINARY
    fprintf(stdout, "Failed reading buffer\n");
  bool binary = false;
  if (!strcmp(buffer, "BINARY")) binary = true;

  if (fscanf(fp, "%s %s", buffer, buffer2) != 2) {
    fclose(fp);
    return 0;
  }

  bool unstructured = false;
  if (!strcmp(buffer, "DATASET") && !strcmp(buffer2, "UNSTRUCTURED_GRID"))
    unstructured = true;

  if ((strcmp(buffer, "DATASET") && strcmp(buffer2, "UNSTRUCTURED_GRID")) ||
      (strcmp(buffer, "DATASET") && strcmp(buffer2, "POLYDATA"))) {
    fprintf(stdout,
            "VTK reader can only read unstructured or poly data datasets\n");
    fclose(fp);
    return 0;
  }

  // read mesh vertices
  int numVertices;
  if (fscanf(fp, "%s %d %s\n", buffer, &numVertices, buffer2) != 3) return 0;
  if (strcmp(buffer, "POINTS") || !numVertices) {
    fprintf(stdout, "No points in dataset\n");
    fclose(fp);
    return 0;
  }
  int data_size;
  if (!strcmp(buffer2, "double"))
    data_size = sizeof(double);
  else if (!strcmp(buffer2, "float"))
    data_size = sizeof(float);
  else {
    fprintf(stdout, "VTK reader only accepts float or double datasets\n");
    fclose(fp);
    return 0;
  }
  fprintf(stdout, "Reading %d points\n", numVertices);
  std::vector<std::array<double, 3>> vertices(numVertices);
  for (int i = 0; i < numVertices; i++) {
    double xyz[3];
    if (binary) {
      if (data_size == sizeof(float)) {
        float f[3];
        if (fread(f, sizeof(float), 3, fp) != 3) {
          fclose(fp);
          return 0;
        }
        if (!bigEndian) SwapBytes((char *) f, sizeof(float), 3);
        for (int j = 0; j < 3; j++) xyz[j] = f[j];
      } else {
        if (fread(xyz, sizeof(double), 3, fp) != 3) {
          fclose(fp);
          return 0;
        }
        if (!bigEndian) SwapBytes((char *) xyz, sizeof(double), 3);
      }
    } else {
      if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
        fclose(fp);
        return 0;
      }
    }
    vertices[i] = {xyz[0], xyz[1], xyz[2]};
    points.push_back(vertices[i]);
  }

  // read mesh elements
  int numElements, totalNumInt;
  if (fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3) {
    fclose(fp);
    return 0;
  }

  bool haveCells = true;
  bool haveLines = false;
  if (!strcmp(buffer, "CELLS") && numElements > 0)
    fprintf(stdout, "Reading %d cells\n", numElements);
  else if (!strcmp(buffer, "POLYGONS") && numElements > 0)
    fprintf(stdout, "Reading %d polygons\n", numElements);
  else if (!strcmp(buffer, "LINES") && numElements > 0) {
    haveCells = false;
    haveLines = true;
    fprintf(stdout, "Reading %d lines\n", numElements);
  } else {
    fprintf(stdout, "No cells or polygons in dataset\n");
    fclose(fp);
    return 0;
  }

  if (haveCells) {
    std::vector<std::vector<uint32_t> > cells(numElements);
    for (auto &cell : cells) {
      int num_vertices, n[100];
      if (binary) {
        if (fread(&num_vertices, sizeof(int), 1, fp) != 1) {
          fclose(fp);
          return 0;
        }
        if (!bigEndian) SwapBytes((char *) &num_vertices, sizeof(int), 1);
        if ((int) fread(n, sizeof(int), num_vertices, fp) != num_vertices) {
          fclose(fp);
          return 0;
        }
        if (!bigEndian) SwapBytes((char *) n, sizeof(int), num_vertices);
      } else {
        if (fscanf(fp, "%d", &num_vertices) != 1) {
          fclose(fp);
          return 0;
        }
        for (int j = 0; j < num_vertices; j++) {
          if (fscanf(fp, "%d", &n[j]) != 1) {
            fclose(fp);
            return 0;
          }
        }
      }
      for (int j = 0; j < num_vertices; j++) {
        if (n[j] >= 0 && n[j] < (int) vertices.size())
          cell.push_back(n[j]);
        else
          fprintf(stdout, "Wrong node index %d\n", n[j]);
      }
    }

    if (unstructured) {
      if (fscanf(fp, "%s %d\n", buffer, &numElements) != 2) {
        fclose(fp);
        return 0;
      }
      if (strcmp(buffer, "CELL_TYPES") || numElements != (int) cells.size()) {
        fprintf(stdout, "No or invalid number of cells types\n");
        fclose(fp);
        return 0;
      }
      for (auto &i : cells) {
        int type;
        if (binary) {
          if (fread(&type, sizeof(int), 1, fp) != 1) {
            fclose(fp);
            return 0;
          }
          if (!bigEndian) SwapBytes((char *) &type, sizeof(int), 1);
        } else {
          if (fscanf(fp, "%d", &type) != 1) {
            fclose(fp);
            return 0;
          }
        }
        switch (type) {
          case 1: elements[0].push_back(i);
            break;//Point
            // first order elements
          case 3: elements[1].push_back(i);
            break;//Line
          case 5: elements[2].push_back(i);
            break;//Triangle
          case 9:elements[3].push_back(i);//Quad
            break;
          case 10:elements[4].push_back(i);//Tet
            break;
          case 12:elements[5].push_back(i);//Hex
            break;
          case 13: elements[6].push_back(i);
            break;//Prism
          case 14: elements[7].push_back(i);
            break;//Pyramid
            // second order elements
          case 21: elements[1].push_back(i);
            break;//Line3
          case 22:elements[2].push_back(i);//Triangle6
            break;
          case 23://quad8
          case 28:elements[3].push_back(i);//quad9
            break;
          case 24:elements[4].push_back(i);//tet10
            break;
          case 25://hex20
          case 29:elements[5].push_back(i);//hex27
            break;
          case 26://Prism15
          case 32:elements[6].push_back(i);//Prism18
            break;
          default: fprintf(stdout, "Unknown type of cell %d\n", type);
            break;
        }
      }
    } else {
      for (auto &cell : cells) {
        int nbNodes = (int) cell.size();
        switch (nbNodes) {
          case 1: elements[0].push_back(cell);
            break;
          case 2: elements[1].push_back(cell);
            break;
          case 3: elements[2].push_back(cell);
            break;
          case 4:elements[3].push_back(cell);
            break;
          default:
            fprintf(stdout,
                    "Unknown type of mesh element with %d nodes\n",
                    nbNodes);
            break;
        }
      }
    }
  } else if (haveLines) {
    if (!binary) {
      int v0, v1;
      char line[100000], *p, *pEnd, *pEnd2;
      for (int k = 0; k < numElements; k++) {

        if (!fgets(line, sizeof(line), fp)) {
          fclose(fp);
          return 0;
        }
        strtol(line, &pEnd, 10);
        v0 = (int) strtol(pEnd, &pEnd2, 10);
        p = pEnd2;
        while (true) {
          v1 = strtol(p, &pEnd, 10);
          if (p == pEnd) break;
          elements[1].push_back({(uint32_t) v0, (uint32_t) v1});
          p = pEnd;
          v0 = v1;
        }
      }
    } else {
      fprintf(stdout, "Line import not done for binary VTK files\n");
    }
  }

  for (auto l : elements[1]) {
    if (l.size() == 2) {
      line.push_back({l[0], l[1]});
    }
  }
  for (auto t : elements[2]) {
    if (t.size() == 3) {
      triangles.push_back({t[0], t[1], t[2]});
    }
  }
  fclose(fp);
  return 1;
}

bool VTKWriterVectorOnPoint(const std::string &meshName,
                            const std::vector<Vec3> &points,
                            const std::vector<ID3> &triangles,
                            std::vector<std::array<double, 9>> &direction_on_triangle) {
  FILE *fp = fopen(meshName.c_str(), "w");
  if (!fp) {
    fprintf(stdout, "Unable to open file '%s'\n", meshName.c_str());
    return false;
  }

  // get the number of vertices and index the vertices in a continuous
  // sequence
  int numVertices = (int) points.size();

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Su\n", SplitFileName(meshName)[1].c_str());
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  // write mesh vertices
  fprintf(fp, "POINTS %d double\n", numVertices);
  for (int i = 0; i < numVertices; i++)
    fprintf(fp, "%.16g %.16g %.16g\n", points[i][0], points[i][1], points[i][2]);

  // loop over all elements we need to save and count vertices

  int numElements = (int) triangles.size(), totalNumInt = (int) triangles.size() * 4;
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);
  for (auto triangle : triangles)
    fprintf(fp, "%d %ld %ld %ld\n", 3, triangle[0], triangle[1], triangle[2]);
  fprintf(fp, "CELL_TYPES %d\n", numElements);
  for (int i = 0; i < numElements; ++i) fprintf(fp, "5\n");

  std::vector<Vec3> direction_on_point(numVertices);

  for (int i = 0; i < triangles.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      direction_on_point[triangles[i][j]] =
          {direction_on_triangle[i][j * 3], direction_on_triangle[i][j * 3 + 1], direction_on_triangle[i][j * 3 + 2]};
    }
  }

  fprintf(fp, "POINT_DATA %d\n", numVertices);
  fprintf(fp, "VECTORS vectors double\n");
  //fprintf(fp, "LOOKUP_TABLE default\n");
  for (int i = 0; i < numVertices; ++i) {
    fprintf(fp, "%.12g %.12g %.12g\n", direction_on_point[i][0], direction_on_point[i][1], direction_on_point[i][2]);
  }

  fclose(fp);

  return true;
}

int main() {
  std::vector<Vec3> points;
  std::vector<ID3> triangles;
  std::vector<ID2> lines;
  std::string infile =
      "D:\\OneDrive\\model\\Mechanical\\vtk_raw_model\\TransmissionHousingFem.vtk";

  readVTK(infile, points, lines, triangles);
  std::vector<Vec3> edge_dir;//边上的标架分量之一
  std::vector<ID3> singularities;
  std::vector<std::array<double, 9>> global_triangle_dir;//三角形的三个点的标架分量之一

  BuildBackgroundMeshAndGuidingField(points,
                                     triangles,
                                     lines,
                                     edge_dir,
                                     global_triangle_dir,
                                     singularities);

  VTKWriterVectorOnPoint("TransmissionHousingFem_vector.vtk",points,triangles,global_triangle_dir);

  return 0;
}