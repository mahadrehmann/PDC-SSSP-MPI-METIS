# Parallel Dynamic SSSP Update Framework using MPI and METIS

🚀 A high-performance parallel framework for updating Single-Source Shortest Paths (SSSP) in large-scale dynamic graphs, based on the paper _"A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks."_ This project is developed as part of our Parallel & Distributed Computing course.

## 📌 Overview

This project explores an efficient two-phase approach to update SSSP trees incrementally in dynamic graphs, instead of recomputing from scratch. We analyze the proposed CPU and GPU implementations, and outline a hybrid parallelization strategy using:

- **MPI** for inter-process communication
- **OpenMP/OpenCL** for intra-node parallelism
- **METIS** for efficient graph partitioning


## 📚 Paper Summary

> **Title:** A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks  
> **Authors:** Xiaojing An, Carl Yang, Ariful Azad, Aydın Buluç  
> **Link:** ([Paper](https://drive.google.com/file/d/1I_W52vP4k3amMnXCPOuUzUgTXi9v3aXj/view))

- Introduces a parallel update framework for dynamic graphs
- Handles batches of edge insertions and deletions
- Updates only affected parts of the graph to save computation
- Implements both shared-memory (OpenMP) and GPU (CUDA) variants
- Proposes VMFB (Vertex-Marking Functional Block) for lock-free parallelism

## 👥 Team Members

- Mahad Rehman Durrani 22i-0792  
- Aniq Noor 22i-0987
- Asim Iqbal 22i-0787

## 🛠️ Tools & Technologies

- MPI (Message Passing Interface)
- OpenMP (shared-memory parallelism)
- METIS (graph partitioning)
- C/C++ for implementation
- PowerPoint / LaTeX for presentation


## 🧠 Phase 1 – Contributions

- Studied and summarized the selected research paper  
- Prepared slides explaining:
  - The core algorithm and update strategy
  - Shared-memory (OpenMP) and GPU (CUDA) implementations
  - Use of METIS for graph partitioning
- Proposed our own implementation strategy using:
  - **MPI + OpenMP**
  - **METIS** for partitioning static graphs before updates

