## An image registration example (a difficult problem with many local minima)
#
# The goal is to discover the translation and rotation needed to align two images.
# You can visualize the images below with
#    using ImageView
#    imshow(img)
#    imshow(imgrot)
# You'll see that they differ by a big rotation.

# You could make this more efficient by using Fourier methods to compute the ideal
# translation given any rotation. But the point here is to illustrate the brute-force approach.

using Images, TestImages, CoordinateTransformations, QuadDIRECT, Base.Test

# Takes a parameter vector `x` and returns a rigid transformation
function tfmx(x, img)
    dx, dy, θ = x
    return Translation(dx, dy) ∘ recenter(RotMatrix(θ), center(img))
end

# Compute the remaining mismatch after applying the transformation implied by `x` to `img`
function mismatch(x, img, target)
    tfm = tfmx(x, img)
    imgw = warp(img, tfm)
    mismatch(imgw, target)
end

function mismatch(img1, img2)
    # Compute the mismatch over the indices that the two images have in common
    inds1, inds2 = indices(img1), indices(img2)
    inds = (intersect(inds1[1], inds2[1]), intersect(inds1[2], inds2[2]))
    (isempty(inds[1]) || isempty(inds[2])) && return Inf
    val = 0.0
    counter = 0
    for x in inds[2], y in inds[1]
	v1, v2 = img1[y, x], img2[y, x]
	if !isnan(v1) && !isnan(v2)
	    val += abs2(v1-v2)
	    counter += 1
	end
    end
    counter == 0 && return Inf
    # This is like "nanmean"
    return val/counter
end

const img = float32.(testimage("lighthouse"))
x = [25, -17, 2*pi/3]
tfm = tfmx(x, img)
const imgrot = warp(img, tfm)

# Test that the ground truth gives a very low mismatch (actually, 0)
@test mismatch(x, img, imgrot) < 1e-3

f(x) = mismatch(x, img, imgrot)

upper = [30, 30, pi] # allow a search of 30 pixels in either direction, plus any rotation
lower = -upper
splits = ([-15, 0, 15], [-15, 0, 15], [-pi/2, 0, pi/2])

# Transformations that differ by less than a pixel or less than 0.02 radians are not
# worth comparing
minwidth = [1, 1, 0.02]

# Simple iterative analysis
root, x0 = analyze(f, splits, lower, upper; maxevals=10, minwidth=minwidth)
# Iterative refinement with control over convergence. `fvalue=0.02` forces it to keep looking until it finds something good.
warn("This will take a long time")
root = analyze!(root, f, x0, splits, lower, upper; maxevals=10000, fvalue=0.02, rtol=0, minwidth=minwidth, print_interval=100)

box = minimum(root)
println("Ground truth: ", x)
println("Solution: ", position(box, x0))
