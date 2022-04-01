
![collapsedvnormal](https://user-images.githubusercontent.com/2405036/161212593-a016f5f9-2be0-4ef6-abca-198fdd9e0814.png)

# Image Driven Simplification

This was the term project for CS453 at Oregon State University which was the implementation of image-driven simplification by Peter Lindstrom and Greg Turk in their paper from 2000.

## Introduction

The goal of this paper is to showcase our implementation of simplification method of a complex 3D image model. The goal of this method is to simplify an image such that it is visually similar, but having the simplified model contain far fewer polygons. Our evaluation of such is done, then, by viewing both the original and simplified model and determining if they are visually similar by our judgement. This is different from other simplification models which use the models edges, vertices, etc. to determine similarity.
We will cover a few previous works that have been done that are similar to this method, then we will cover our methods and division of tasks for simplifying the Stanford bunny model. Finally, we will cover the evaluation and interpretation of our simplified model and determine success or failure of such.

## Background

There is a lot of benefit and reasoning for simplifying a highly complex 3-D model. Firstly, in the context of rendering for industries such as video game or applications such as maps, decreasing complexity logically decreases load and processing time, which is important in the context of video games or visual applications. Time to load and storage are important when there is either limited processing or storage; such as on a phone.
It is worth noting, there are applications that need geometric exactness, which would not be applicable to this simplification method. However, there are many simulations that only require visual similarity, those include "vehicle simulations, building walk-through, video games, acceleration of off-line rendering, and fast manipulation of simplified models that can then be rendered at high fidelity after the user has chosen a suitable view". In many aspects of image-rendering and usage, simplifying a model while maintaining it's visual representation can increase speed of load and visualization, decrease size of the model for storage, and cause less load in processing.

## Goals

- Render and load 3D mesh data of the Stanford Bunny
- Implement edge-cost function to determine the cost of removing each edge
- Implement edge-remove function that collapses *n* of the lowest-costing edges
- Render the model and compare it to the original to determine if it can be simplified further or not

![fullcollapse_bunny](https://user-images.githubusercontent.com/2405036/161212657-66efa3ce-4959-4d4b-a2c1-301d95c04b3c.png)

## Results

Once the program was completed, the Stanford Bunny Mesh was simplified by a specified n number, which is 1511 triangles. On average, we collapse around 15% of the total number of triangles in the model. The edge collapse algorithm worked properly and detected which edges had the lowest impact, and destroyed the edge. Additionally, any edges that were detected to be on a boundary, which would have a high visual impact if collapsed, were always ignored by the collapse algorithm to avoid a great change visually in the model. The incomplete part of the program, however, is that the rendered image did not have reconstructed edges, so there are gaps in the model where an edge was collapsed. With more time, the reconstruction of the Stanford bunny would be able to be completed. However, since this project was completed within the constraints of an assignment with a specified due date, it was not fully completed.

## Evaluation and Interpretation

As mentioned above, comparing the original model to the simplified model, it can be seen that the gaps where an edge collapsed was in a low-impact area where, if reconstructed, would have very minimal visual difference compared to the original model. Seen in figure 2, the ~3 collapsed edges below the bunny's chin, the chest area is fairly flat, so a collapse and smoothing of that area would not be greatly noticeable to the human eye. With just collapsing of different triangles which are deemed not worth in the terms of the algorithm can lead to faster load times in the future with proper reconstructed edges. It takes our algorithm 95 iterations to collapse all 1511 triangles which is extremely fast. As you can see in the full edge collapse in figure 2, we were able to eliminate unnecessary vertices in the smoother parts of the model. Significant amount of resources and memory is required though for this algorithm to be able to compute due to the amount of matrices and vectors that need to be computed to determine edge cost. However, we face another issue where on first run the program crashes due to random initialization of the data set leading the user to think the project is incomplete, this is not the case and the project can just be re-ran without issue.

## Conclusion

In conclusion, our goal for the project, to correctly implement and edge collapse algorithm to simplify an image was successful. To determine how many edges to collapse, the cost per each edge is calculated and through applying certain constraints, edges are collapsed in order of lowest cost. The final rendering of the model with simplified and re-constructed edges was not completed, but we hope to be able to apply more concepts from Lindstrom and Turk's paper to be able to fully render a simplified Stanford Bunny model.

## References

[1] P. Lindstrom and G. Turk, “Image-driven simplification,” vol. 19, no. 3,
pp. 204–241.

[2] W. J. Schroeder, J. A. Zarge, and W. E. Lorensen, “Decimation of triangle
meshes,” Proceedings of the 19th annual conference on Computer graphics
and interactive techniques - SIGGRAPH 92, p. 65–70, Jul 1992.

[3] R. Ronfard and J. Rossignac, “Full-range approximation of triangulated
polyhedra.,” Computer Graphics Forum, vol. 15, p. 67–76, Aug 1996.

[4] S.-e. Yoon, B. Salomon, M. Lin, and D. Manocha, “Fast collision detection
between massive models using dynamic simplification.,” pp. 139–150, 01
2004.

[5] X. Wu and J. Ye, “A new 3-d mesh simplificati  ́on algorithm.,” Revista
Colombiana de Computaci ́on, vol. 4, 01 2003.

[6] K. Ito, “Mesh simplification in scaniverse,” Jun 2021.
