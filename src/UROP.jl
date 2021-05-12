using Images, ImageView, FileIO, StatsBase, Statistics, CSV, Tables, Glob, LinearAlgebra
import Base.Threads

function ReadTIFF(filename,x,y)
    fid = open(filename, "r")

    T = Array{UInt16}(undef, x, y)
    read!(fid,T)
    close(fid);
    T = convert(Array{Float32}, T)./(2^16-1)
    return T
end

function ReadStack(dir, filter=nothing)
    myDir = readdir(glob"*rec000*.tif",dir)
    size = Images.size(Images.load(myDir[1]))
    N = length(myDir)
    stack = Array{Float32}(undef, size[2],size[1], N)
    Threads.@threads for i in 1:N
        stack[:,:,i]=ReadTIFF(myDir[i],size[2],size[1])
        if filter != nothing
            stack[:,:,i] = imfilter(stack[:,:,i], filter)
        end
    end
    return stack
end

function reformat(array)
    x, y, z = size(array)
    result = []
    for i in range(1,length=x)
        for j in range(1,length=y)
            for k in range(1,length=z)
                if array[i,j,k] != 0
                    push!(result, array[i,j,k])
                end
            end
        end
    end
    return result
end

function find_ribs(bones, lookahead = 50, lookback = 30)
    ribs = BitArray{3}(undef,size(bones)...)
    N = size(bones, 3)
    total = sum(view(bones,:,:,1:min(lookahead,N)),dims=3)
    ribs[:,:,1:1] .= (0 .< total)
    for k in 2:N
        if k-lookback >= 2
            total .-= view(bones, :, :, k-lookback-1)
        end
        if k+lookahead <= N
            total .+= view(bones, :, :, k+lookahead)
        end
        ribs[:,:,k:k] .= (0 .< total)
    end
    ribs = label_components(ribs);
    ribs = ribs.==find_biggest_component(ribs);
    return ribs
end

function find_biggest_component(y)
    u=unique(y)
    d=Dict([(i,count(x->x==i,y)) for i in u])
    max = 1
    for key in keys(d)
        if key != 0 && d[key]>d[max]
            max = key
        end
    end
    return max
end

function find_centroids(ribs)
    cm = zeros((552,472,530))
    data = Dict{Int,Array{Int64,1}}()
    for k in range(20,length=390)
        mx= 0.0
        my=0.0
        m = 0.0
        for j in range(1,length=472)
            for i in range(1,length=552)
                 mx += ribs[i,j,k] * i;
                 my += ribs[i,j,k] * j;
                 m += ribs[i,j,k];
            end
        end
        cmx = Int(round(mx/m))
        cmy = Int(round(my/m))
        cm[cmx,cmy,k] = 1
        get!(data,k,[cmx,cmy])
    end
    return cm,data
end

function apply_mask(mask, img)
    result_img = RGB.(Gray.(img))
    (width, height, depth) = size(img)
    for k in 1:depth
        for j in 1:height
            for i in 1:width
                if mask[i,j,k] == 1
                    result_img[i,j,k] = ARGB(RGB(1, 0, 0),0.01)
                end
            end
        end
    end
    return result_img
end

function grid_points_in_poly(shape, verts)
    @assert length(shape) == 2
    out = Array{Bool}(undef, shape)
    for I in CartesianIndices(out)
        @inbounds out[I] = point_in_polygon(verts, I.I...)
    end
    out
end

function point_in_polygon(verts::Vector{CartesianIndex{2}}, x, y, eps=1e-12)
    # Initialize the loop
    p1 = verts[end]
    x1, y1 = p1.I
    x1 -= x
    y1 -= y
    r_cross, l_cross = UInt(0), UInt(0)
    for point in verts
        xp, yp = point.I
        x0 = xp - x
        y0 = yp - y
        if (-eps < x0 < eps) && (-eps < y0 < eps)
            # it is a vertex with an eps tolerance
            return true # VERTEX
        end

        # if e straddles the x-axis
        if ((y0 > 0) != (y1 > 0))
            # check if it crosses the ray
            if ((x0 * y1 - x1 * y0) / (y1 - y0)) > 0
                r_cross += 1
            end
        end
        # if reversed e straddles the x-axis
        if ((y0 < 0) != (y1 < 0))
            # check if it crosses the ray
            if ((x0 * y1 - x1 * y0) / (y1 - y0)) < 0
                l_cross += 1
            end
        end
        x1 = x0
        y1 = y0
    end
    if (r_cross & 1) != (l_cross & 1)
        # on edge if left and right crossings not of same parity
        return true # EDGE
    end
    if (r_cross & 1) == 1
        # inside if odd number of crossings
        return true # INSIDE
    end
    # outside if even number of crossings
    return false # OUTSIDE
end

function ROI(ribs)
    y = ImageMorphology.label_components(ribs);
    c = find_biggest_component(y);
    lungs = y.== c
    Threads.@threads for i in range(1,length=size(lungs)[3])
        chull = convexhull(lungs[:,:,i])
        lungs[:,:,i] = grid_points_in_poly(size(lungs[:,:,i]), chull)
    end
    return lungs
end


function calculate_data(dir)
    if occursin(r"_Rec",dir)
        mydir=[dir]
    else
        cd(dir)
        mydir=glob("*_Rec")
    end
    result = ["Name" "Mean lungs" "Median Lungs" "Mean bones" "Median bones" "Mean xgrad lungs" "Median xgrad lungs"  "Mean ygrad lungs" "Median ygrad lungs"  "Mean zgrad lungs" "Median zgrad lungs";]
    # result = ["Name" "Mean lungs" "Median Lungs" "Mean bones" "Median bones" ];
    ch = Channel{Any}(Inf)
    try
        for i in range(1, length=length(mydir))
            name = joinpath(dir,mydir[i])
            mouse = ReadStack(name)
            bones = mouse .> .68
            ribs = find_ribs(bones)
            in_ribs = ROI(ribs)
            # lungs
            interior = mouse.*in_ribs
            mod_interior = reformat(interior)
            # bones
            bones = reformat(bones.*mouse)
            # # gradients
            x =gradx(interior)
            y =grady(interior)
            z =gradz(interior)
            gx = reformat(x)
            gy = reformat(y)
            gz = reformat(z)
            #here is where we would add all other data calculations
            row = [mydir[i] string(Statistics.mean(mod_interior)) string(Statistics.median(mod_interior)) string(Statistics.mean(bones)) string(Statistics.median(bones)) string(Statistics.mean(gx)) string(Statistics.median(gx))  string(Statistics.mean(gy)) string(Statistics.median(gy))  string(Statistics.mean(gz)) string(Statistics.median(gz))]
            # row = [mydir[i] string(Statistics.mean(mod_interior)) string(Statistics.median(mod_interior)) string(Statistics.mean(bones)) string(Statistics.median(bones))]
            put!(ch, row)

        end
    catch
        close(ch)
        result = reduce(vcat, ch,init=result)
        t=Tables.table(result)
        write_to = joinpath(dir,"Data.csv")
        CSV.write(write_to,t)
        return
    end
    close(ch)
    result = reduce(vcat, ch,init=result)
    t=Tables.table(result)
    write_to = joinpath(dir,"Data.csv")
    CSV.write(write_to,t)
end

function gradx(array)
    N = size(array,1)
    grad = Array{Float32}(undef,size(array)...)
    for i in range(1,length=N)
        if i != 1 && i!= N
            @views grad[i,:,:] =(array[i+1,:,:] .- array[i-1,:,:])./2
        elseif i== 1
            @views grad[i,:,:] =(array[i+1,:,:] .- array[i,:,:])./2
        else
            @views grad[i,:,:] =(array[i,:,:] .- array[i-1,:,:])./2
        end
    end

    return grad
end

function grady(array)
    N = size(array,2)
    grad = Array{Float32}(undef,size(array)...)
    for i in range(1,length=N)
        if i != 1 && i!= N
            @views grad[:,i,:] =(array[:,i+1,:] - array[:,i-1,:])./2
        elseif i== 1
            @views grad[:,i,:] =(array[:,i+1,:] - array[:,i,:])./2
        else
            @views grad[:,i,:] =(array[:,i,:] - array[:,i-1,:])./2
        end
    end

    return grad
end

function gradz(array)
    N = size(array,3)
    grad = Array{Float32}(undef,size(array)...)
    for i in range(1,length=N)
        if i != 1 && i!= N
            @views grad[:,:,i] =(array[:,:,i+1] - array[:,:,i-1])./2
        elseif i== 1
            @views grad[:,:,i] =(array[:,:,i+1] - array[:,:,i])./2
        else
            @views grad[:,:,i] =(array[:,:,i] - array[:,:,i-1])./2
        end
    end

    return grad
end

function save_image(roi,dir)
    N = size(roi,3)
    new_dir = joinpath(dir,"images")
    for i in range(1,length=N)
        new_file = basename(dir)*"_layer$i.png"
        save(joinpath(new_dir,new_file), colorview(Gray,(roi[:,:,i])))
    end
end

function calculate_singular_hist(mouse, region_of_interest=false,save_images=false,dir=false)
    if region_of_interest
        bones = mouse.>.68
        ribs = find_ribs(bones)
        roi = mouse.*ROI(ribs)
    else
        roi = mouse
    end
    if save_images
        save_image(roi,dir)
    end

    h = fit(StatsBase.Histogram, vec(roi),nbins=5000)
    result = h.edges[1][3:length(h.edges[1])-1]
    result = hcat(result,h.weights[2:end-1])

    if region_of_interest
        new_file = basename(dir)*"_roi_histogram.csv"
    else
        new_file = basename(dir)*"_histogram.csv"
    end

    CSV.write(joinpath(dir,new_file),Tables.table(result))
    return result
end

function preprocess_data(dir)
    mouse = ReadStack(dir)
    bones = mouse.>.68
    ribs = find_ribs(bones)
    roi = mouse.*ROI(ribs)
    return roi
end


function calculate_histograms(dir, region_of_interest=false, save_images=false)
    print("Calculating many histograms")
    cd(dir)
    mydir=glob("*_Rec")
    result = ["Name"]
    ch = Channel{Any}(Inf)
    Threads.@threads for i in 1:length(mydir)
        print(i)
        name = joinpath(dir,mydir[i])
        if region_of_interest
            mouse = preprocess_data(name)
        else
            mouse = ReadStack(name)
        end

        t=calculate_singular_hist(mouse, region_of_interest,save_images,name)
        if i==1
            result = vcat(result, t[:,1])
        end
        column = [mydir[i]]
        column = vcat(column,t[:,2])
        put!(ch, column)

    end
    close(ch)
    result = reduce(hcat, ch,init=result)

    t=Tables.table(result)
    if region_of_interest
        new_file = basename(dir)*"_roi_histograms.csv"
    else
        new_file = basename(dir)*"_histograms.csv"
    end

    CSV.write(joinpath(dir,new_file),t)

    return
end
