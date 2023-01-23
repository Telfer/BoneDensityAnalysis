library(tidyverse)
library(rdist)
library(oro.nifti)


#' Import Surface PLY file
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param filepath String.
#' @param ply_type String. Which software generated the mesh. Currently, "Meshlab"
#'   and "Rhino" supported
#' @return List. List containing vertex coordinates, vertex normals, vertex
#'   colors, and face data
import_ply <- function(filepath, readColor = FALSE, ply_type = "Meshlab") {
  # load file as text strings
  x <- readLines(filepath)
  
  # determine number of vertices
  vertex_n_line <- str_which(x, "element vertex")
  vertex_n <- as.numeric(str_extract(x[vertex_n_line], "(\\d)+"))
  
  # header
  start_line <- str_which(x, "end_header")
  header <- x[1:start_line]
  
  # import ply vertex data as data frame
  vertex_data <- x[(start_line + 1):(start_line + vertex_n)]
  tc_vertex_data <- textConnection(vertex_data)
  vertex_data <- read.table(tc_vertex_data, sep = " ")
  close(tc_vertex_data)
  vertex_coords <- vertex_data[, 1:3]
  vertex_normals <- vertex_data[, 4:6]
  if (readColor == TRUE) {vertex_colors <- vertex_data[, 7:9]}
  
  # import ply face data as data frame
  face_data <- x[(start_line + vertex_n + 1):length(x)]
  tc_face_data <- textConnection(face_data)
  face_data <- read.table(tc_face_data, sep = " ")
  close(tc_face_data)
  face_data <- face_data[, 1:4] # add check here to check mesh is triangular
  
  # return
  if (readColor == TRUE) {
    output <- list(header = header, vertex_coords = vertex_coords, 
                   vertex_normals = vertex_normals, face_data = face_data, 
                   vertex_colors = vertex_colors)
  } else {
    output <- list(header = header, vertex_coords = vertex_coords, 
                   vertex_normals = vertex_normals, face_data = face_data)
  }
  
  # return
  return(output)
}


#' check landmarks are close to the bone
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh_path String. Filepath to triangulated surface mesh in 
#' ascii ply format
#' @param landmark_path String. Filepath to landmark data in .fcsv format from 
#' 3D Slicer 
#' @param threshold Numeric. Distance landmark can be from surface without 
#' warning being thrown
#' @return String. Returns a message warning that landmarks are not on bone 
#' surface
landmark_check <- function(surface_mesh_path, landmark_path, threshold = 1.0) {
  surface_mesh <- import_ply(surface_mesh_path)
  vertices <- surface_mesh$vertex_coords
  landmarks <- read.csv(landmark_path, skip = 2)[2:4]
  
  dists <- c()
  for (i in 1:nrow(landmarks)) {
    x <- cdist(vertices, landmarks[i, ])
    dists <- c(dists, min(x))
  }
  
  # return message if landmarks not on bone surface
  if (any(dists > threshold)) {print("landmarks not on bone surface")}
}


#' Sigma beta CT calculations
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param table_height
#' @param calibration_curves
#' @param scanner "CT1" or "CT2"
#' @return Vector with two elements, sigma and beta
ct_coefficients <- function(table_height, calibration_curves, scanner, return_coeff = "sigma") {
  # get curves for scanner
  calibration_curves <- calibration_curves %>% filter(Scanner == scanner)
  sigmaCT <- approx(calibration_curves$TableHeight, calibration_curves$sigma, table_height)$y
  betaCT <- approx(calibration_curves$TableHeight, calibration_curves$beta, table_height)$y
  
  # return
  if(return_coeff == "sigma") {return(sigmaCT)}
  if(return_coeff == "beta") {return(betaCT)}
}


#' Redefine surface points
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adapted from geomorph
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param landmarks Data frame. Contains 3D coords of landmarks
#' @param no_surface_sliders Numeric. Number of surface points to add 
#' @return Data frame. 3D coords of remapped surface points
surface_points_template <- function(surface_mesh, landmarks, no_surface_sliders) {
  # get points from surface
  vertices <- as.matrix(surface_mesh$vertex_coords)
  colnames(vertices) <- c("xpts", "ypts", "zpts")
  
  # landmarks
  lmk.add <- NULL
  for(i in 1:nrow(landmarks)){
    lmk.add <- rbind(lmk.add, which.min(sqrt((landmarks[i, 1] - vertices[, 1]) ^ 2 + 
                                             (landmarks[i, 2] - vertices[, 2]) ^ 2 + 
                                             (landmarks[i, 3] - vertices[, 3]) ^ 2))[1])}
  nlandmarks <- nrow(landmarks) 
  vertices <- vertices[-lmk.add, ]
  
  # calculate new points
  colnames(landmarks) <- c("xpts", "ypts", "zpts")
  new_surface_points <- rbind(landmarks, 
                              kmeans(x = vertices, centers = no_surface_sliders, 
                                     iter.max = 25)$centers)
  
  # return
  return(new_surface_points)
}


#' New surface points from template
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adpated from geomorph
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param landmarks Data frame. Contains 3D coords of landmarks 
#' @param template Data frame. 3D coords of remapped surface points
#' @return Data frame. 3D coords of remapped surface points
surface_points_new <- function(surface_mesh, landmarks, template) {
  ## helper functions
  rotate.mat <- function(M, Y){
    k <- ncol(M)
    M <- cs.scale(M); Y <- cs.scale(Y)
    MY <- crossprod(M, Y)
    sv <- La.svd(MY, k, k)
    u <- sv$u; u[, k] <- u[, k] * determinant(MY)$sign
    v <- t(sv$vt)
    tcrossprod(v, u)
  }
  
  csize <- function(x) sqrt(sum(center(as.matrix(x))^2))
  
  center <- function(x){
    if(is.vector(x)) x - mean(x) else {
      x <- as.matrix(x)
      dims <- dim(x)
      fast.center(x, dims[1], dims[2])
    }
  }
  
  fast.center <- function(x, n, p){
    m <- colMeans(x)
    x - rep.int(m, rep_len(n, p))
  }
  
  cs.scale <- function(x) x/csize(x)
  
  fast.solve <- function(x) { 
    x <- as.matrix(x)
    if(det(x) > 1e-8) {
      res <- try(chol2inv(chol(x)), silent = TRUE)
      if(class(res) == "try-error") res <- fast.ginv(x)
    } else res <- fast.ginv(x)
    return(res)
  }
  
  tps2d3d <- function(M, matr, matt, PB = TRUE){		#DCA: altered from J. Claude 2008
    p <- dim(matr)[1]; k <- dim(matr)[2]; q <- dim(M)[1]
    Pdist <- as.matrix(dist(matr))
    ifelse(k == 2, P <- Pdist^2*log(Pdist^2), P <- Pdist)
    P[which(is.na(P))] <- 0
    Q <- cbind(1, matr)
    L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0, k + 1, k + 1)))
    m2 <- rbind(matt, matrix(0, k + 1, k))
    coefx <- fast.solve(L)%*%m2[, 1]
    coefy <- fast.solve(L)%*%m2[, 2]
    if(k == 3){coefz <- fast.solve(L)%*%m2[, 3]}
    fx <- function(matr, M, coef, step){
      Xn <- numeric(q)
      for (i in 1:q){
        Z <- apply((matr-matrix(M[i,], p, k, byrow = TRUE))^2, 1, sum)
        ifelse(k == 2, Z1<-Z*log(Z), Z1<-sqrt(Z)); Z1[which(is.na(Z1))] <- 0
        ifelse(k == 2, Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + sum(coef[1:p]*Z1),
               Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + coef[p+4]*M[i,3] + sum(coef[1:p]*Z1))
        if(PB == TRUE){setTxtProgressBar(pb, step + i)}
      }
      return(Xn)
    }
    matg <- matrix(NA, q, k)
    if(PB==TRUE){pb <- txtProgressBar(min = 0, max = q*k, style = 3) }
    matg[,1] <- fx(matr, M, coefx, step = 1)
    matg[,2] <- fx(matr, M, coefy, step=q)
    if(k==3){matg[,3] <- fx(matr, M, coefz, step=q*2)
    }
    if(PB==TRUE) close(pb)
    return(matg)
  }
  
  fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
    X <- as.matrix(X)
    k <- ncol(X)
    Xsvd <- La.svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <-((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <-t(Xsvd$vt)[, Positive, drop = FALSE]
    v%*%rtu
  }
  
  # format
  bone <- as.matrix(surface_mesh$vertex_coords)
  
  # closet vertex to landmark
  lmk.add <- NULL
  for(i in 1:nrow(landmarks)){
    lmk.add <- rbind(lmk.add, 
                     which.min(sqrt((landmarks[i, 1] - bone[, 1]) ^ 2 + 
                                    (landmarks[i, 2] - bone[, 2]) ^ 2 + 
                                    (landmarks[i, 3] - bone[, 3]) ^ 2))[1])
  } 
  
  nlandmarks <- nrow(landmarks)
  
  # center bone
  bone_centered <- center(bone)
  bone_trans <- colMeans(bone)
  
  # center template
  template <- center(template) * (csize(bone_centered[lmk.add, ]) / csize(template[(1:nlandmarks), ]))  
  template <- template %*% rotate.mat(bone_centered[lmk.add, ], template[(1:nlandmarks), ])
  
  # sliding points
  template.tps <- tps2d3d(template[-(1:nlandmarks), ], template[(1:nlandmarks), ], bone_centered[lmk.add, ])             
  spec.surfs <- bone_centered[-lmk.add, ]
  nei <- numeric(dim(template.tps)[1])
  sliders <- matrix(NA, nrow = dim(template.tps)[1], ncol = 3)
  for (i in 1:dim(template.tps)[1])     {
    nei[i] <- which.min(sqrt((template.tps[i, 1] - spec.surfs[, 1]) ^ 2 + 
                             (template.tps[i, 2] - spec.surfs[ ,2]) ^ 2 + 
                             (template.tps[i, 3] - spec.surfs[, 3]) ^ 2))[1]
    sliders[i,] <- spec.surfs[nei[i], ]
    spec.surfs <- spec.surfs[-nei[i], ]  
  }
  
  # make output matrix
  selected.out <- rbind(bone_centered[lmk.add, ], sliders)
  
  # translate back
  selected.out[, 1] <- selected.out[, 1] - (bone_trans[1] * - 1)
  selected.out[, 2] <- selected.out[, 2] - (bone_trans[2] * - 1)
  selected.out[, 3] <- selected.out[, 3] - (bone_trans[3] * - 1)
  
  # return
  return(selected.out)
}


#' Finds material properties of bone at point
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param mapped_coords Data frame. 3D coords of remapped surface points
#' @param normal_dist Numeric. Distance surface normal should penetrate surface
#' @param nifti Nifti
#' @param betaCT Numeric. Calibration value for CT to density calculation
#' @param sigmaCT. Numeric. Calibration value for CT to density calculation
#' @param voxel_size Numeric. Defaults to 0.5. Currently only isotropic spacing
#' @return Vector. Vector with value for each point on surface
surface_normal_intersect <- function(surface_mesh, mapped_coords, normal_dist = 3.0, nifti, 
                                     betaCT = 1.0, sigmaCT = 1.0, 
                                     voxel_size = 0.5) {
  # format surface data
  surface_coords <- as.matrix(surface_mesh$vertex_coords)
  surface_normals <- as.matrix(surface_mesh$vertex_normals)
  
  # format new point data
  vertex_coords <- data.matrix(mapped_coords)
  dims <- dim(vertex_coords)
  vertex_coords <- as.numeric(vertex_coords)
  dim(vertex_coords) <- dims 
  
  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_seq <- rev(seq(slot(nifti, "qoffset_x"), by = voxel_size, length.out = dims[1]))
  y_seq <- seq(slot(nifti, "qoffset_y"), by = voxel_size, length.out = dims[2])
  z_seq <- seq(slot(nifti, "qoffset_z"), by = voxel_size, length.out = dims[3])
  
  # check bone is within scan volume
  bone_x_min <- min(vertex_coords[, 1])
  bone_x_max <- max(vertex_coords[, 1])
  bone_y_min <- min(vertex_coords[, 2])
  bone_y_max <- max(vertex_coords[, 2])
  bone_z_min <- min(vertex_coords[, 3])
  bone_z_max <- max(vertex_coords[, 3])
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)
  x1_good <- bone_x_min > vol_x_min
  x2_good <- bone_x_max < vol_x_max
  y1_good <- bone_y_min > vol_y_min
  y2_good <- bone_y_max < vol_y_max
  z1_good <- bone_z_min > vol_z_min
  z2_good <- bone_z_max < vol_z_max
  vals <- c(x1_good, x2_good, y1_good, y2_good, z1_good, z2_good)
  print(c("x: ", bone_x_min, bone_x_max, vol_x_min, vol_x_max))
  print(c("y: ", bone_y_min, bone_y_max, vol_y_min, vol_y_max))
  print(c("z: ", bone_z_min, bone_z_max, vol_z_min, vol_z_max))
  if (all(vals) != TRUE) {stop("bone not within scan volume")}
  
  # Find voxels intercepted by line
  mat_peak <- rep(NA, times = nrow(vertex_coords))
  for (i in 1:nrow(vertex_coords)) {
    # find nearest point
    yy <- t(as.matrix(vertex_coords[i, ], 1, 3))
    y <- cdist(yy, surface_coords)
    matched_point <- which.min(y)
    
    # points to test
    start_point <- surface_coords[matched_point, ]
    end_point <- surface_coords[matched_point, ] + (surface_normals[matched_point, ] * -1 * normal_dist)
    px <- seq(from = start_point[1], to = end_point[1], length.out = 10)
    py <- seq(from = start_point[2], to = end_point[2], length.out = 10)
    pz <- seq(from = start_point[3], to = end_point[3], length.out = 10)
    
    # for each point along line find voxel with max value
    max_line <- rep(NA, times = 10)
    for (j in 1:10) {
      voxel <- c(which.min(abs(px[j] - x_seq)),
                 which.min(abs(py[j] - y_seq)),
                 which.min(abs(pz[j] - z_seq)))
      
      max_line[j] <- img_data[voxel[1], voxel[2], voxel[3]]
    }
    
    # add HU column
    mat_peak[i] <- max(max_line)
  }
  
  # convert to density
  mat_peak <- ((mat_peak - betaCT) / sigmaCT)
  
  # return
  return(mat_peak)
}


#' Finds point closest to vertex for all vertices in a surface mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertex_coords Matrix. Generated from import_ply function
#' @param template_points matrix
#' @return Vector. Closest point on mesh to each template pint
mesh_template_match <- function(surface_mesh, template_points) {
  # get vertex coords
  vertex_coords <- as.matrix(surface_mesh$vertex_coords[,-4])
  
  # calculate distances
  y <- cdist(vertex_coords, template_points)
  
  # calculate which point is closest
  matched_points <- apply(y, 1, which.min)
  
  # return
  return(matched_points)
}


#surface_mesh_path <- "Z:/Projects/058_Pelvis Bone Quality/data/PBQ002/surface_models/pelvisLeft_surfaceModel1-5.ply"
#landmark_path <- "Z:/Projects/058_Pelvis Bone Quality/data/PBQ002/surface_models/PBQ002_landmarks_right.fcsv"
#surface_mesh <- import_ply(surface_mesh_path)
#landmarks <- read.csv(landmark_path, skip = 2)[, 2:4]
