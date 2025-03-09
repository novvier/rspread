subroutine new_barrier(xyzsrc, nrow, ncol, size, dem, land, off, roff, barheight, bardist, vegmax)
    implicit none

    integer(kind=4), parameter :: npts = 100
    integer(kind=4), intent(in) :: nrow, ncol
    integer(kind=4), intent(in) :: land(nrow,ncol)
    integer(kind=4) :: m, n, i
    integer(kind=4) :: lc(npts)
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), parameter :: m2ft = 3.28084
    real(kind=r8), intent(in) :: xyzsrc(3), size, dem(nrow,ncol)
    real(kind=r8), dimension(nrow,ncol), intent(out) :: barheight, bardist, vegmax
    real(kind=r8), intent(in) :: off, roff ! offset, roffset
    real(kind=r8) :: max_dist, src_elev, xyzrec(3)
    real(kind=r8) :: hgt(npts), dist_vec(npts)
    real(kind=r8) :: distance, slope, bar_dist, slope_elev, max_veg_value
    real(kind=r8) :: ha_slope, max_height, receiver_elev

    ! PARAMETROS DE ENTRADA
    !size = 30.48
    !xyzsrc(1) = 58.0
    !xyzsrc(2) = 47.0
    !roffset = 1
    !offset = 0.34
    !xyzsrc(3) = dem(int(xyzsrc(2)), int(xyzsrc(1)))

    src_elev = xyzsrc(3) + off * m2ft

    ! Evaluar para cada celda
    do m = 1, ncol
        do n = 1, nrow
            ! Definir las coordenadas de la celda
            xyzrec(1) = m
            xyzrec(2) = n
            xyzrec(3) = dem(n, m)

            ! Llamar a la subrutina para obtener las coordenadas del corte del terreno
            call get_terrain_cut(xyzsrc, xyzrec, size, npts, ncol, nrow, dem, land, max_dist, dist_vec, hgt, lc)

            distance = max_dist * m2ft
            if (distance == 0) then
                distance = 0.01
            end if

            receiver_elev = xyzrec(3) + roff * m2ft

            slope = (receiver_elev - src_elev) / distance

            max_height = 0
            bar_dist = distance

            do i = 1, npts
                slope_elev = slope * dist_vec(i) * m2ft + src_elev
                ha_slope = hgt(i) - slope_elev
                if (ha_slope > max_height) then
                    max_height = ha_slope
                    bar_dist = m2ft * dist_vec(i)
                end if
            end do

            barheight(n,m) = max_height
            bardist(n,m) = bar_dist

            call calc_vegmax(lc, distance, npts, m2ft, max_veg_value)
            vegmax(n,m) = max_veg_value
        end do
    end do

    print *, 'Results iterative calculation successful'

end subroutine new_barrier

subroutine get_terrain_cut(xysrc, xyrec, size, npts, ncol, nrow, mdem, mland, max_dist, dist_vec, hgt, lc)
    implicit none

    integer(kind=4), intent(in) :: npts, ncol, nrow
    integer(kind=4), intent(in) :: mland(nrow,ncol)
    integer(kind=4), intent(out) :: lc(npts)
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: size, mdem(nrow,ncol), xysrc(2), xyrec(2)
    real(kind=r8), intent(out) :: max_dist, dist_vec(npts), hgt(npts)
    real(kind=r8) :: rland(nrow, ncol), rlc(npts)
    real(kind=r8) :: x_coords(npts) , y_coords(npts)
    real(kind=r8) :: xdiff, ydiff, xsource, ysource, xreceiver, yreceiver, minvalue

    xsource = xysrc(1)
    ysource = xysrc(2)
    xreceiver = xyrec(1)
    yreceiver = xyrec(2)

    xdiff = (xreceiver - xsource) * size
    ydiff = (yreceiver - ysource) * size
    max_dist = sqrt(xdiff**2 + ydiff**2)

    call find_distances(max_dist, npts, dist_vec)

    ! Get x_coords
    call find_coords(xsource, xreceiver, npts, x_coords)
    ! Get y_coords
    call find_coords(ysource, yreceiver, npts, y_coords)

    ! Comprobación de valores mínimos:
    call min(x_coords, minvalue, npts)

    if(floor(minvalue) < 1) then
        print *, 'Error: x_coords is less than 1 for', xreceiver, yreceiver
        print *, 'x_coords:'
        print *, x_coords
        error stop
    end if

    call min(y_coords, minvalue, npts)

    if(floor(minvalue) < 1) then
        print *, 'Error: y_coords is less than 1 for', xreceiver, yreceiver
        print *, 'y_coords:'
        print *, y_coords
        error stop
    end if

    ! use coordinates to extract elevation
    call extract_values(ncol, nrow, mdem, npts, x_coords, y_coords, 'bilinear', hgt)

    ! use coordinates to extract land cover
    rland = real(mland, kind=8)
    call extract_values(ncol, nrow, rland, npts, x_coords, y_coords, 'nearests', rlc)
    lc = nint(rlc)

end subroutine get_terrain_cut

subroutine find_distances(max_dist, npts, dist_vec)
    implicit none

    integer(kind=4), intent(in) :: npts
    integer(kind=4) :: i
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: max_dist
    real(kind=r8), intent(out) :: dist_vec(npts)
    real(kind=r8) :: increment

    increment = max_dist / real(npts - 1, kind=8)

    if (increment == 0.0) then
        dist_vec = 1.0
    else

        do i = 1, npts
            dist_vec(i) = (i - 1) * increment
            dist_vec(npts) = max_dist
        end do
    end if

end subroutine find_distances

subroutine find_coords(source, receiver, npts, coords_1d)
    implicit none

    integer(kind=4), intent(in) :: npts
    integer(kind=4) :: i
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: source, receiver
    real(kind=r8), dimension(npts), intent(out)  :: coords_1d
    real(kind=r8), dimension(npts) :: coords_ini
    real(kind=r8) :: abs_diff, direction, diff

    diff = receiver - source

    direction = sign(real(1.0, kind=8), diff)
    abs_diff = abs(diff)

    call find_distances(abs_diff, npts, coords_ini)

    do i = 1, npts
        coords_1d(i) = source + coords_ini(i) * direction
    end do

end subroutine find_coords

subroutine extract_values(ncol, nrow, in_array, npts, x_coords, y_coords, extract_type, values)
    implicit none

    integer(kind=4), intent(in) :: ncol, nrow, npts
    integer(kind=4) :: i, y, x
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: in_array(nrow, ncol), x_coords(npts), y_coords(npts)
    real(kind=r8), intent(out) :: values(npts)

    character (len=8), intent(in) :: extract_type

    do i = 1, size(x_coords)
        if (extract_type == 'nearests') then
            x =  floor(x_coords(i))
            y =  floor(y_coords(i))
            values(i) = in_array(y, x)
        else if (extract_type == 'bilinear') then
            ! Bilinear interpolation
            call do_bilinear_interpolation(ncol, nrow, in_array, x_coords(i), y_coords(i), values(i))
        end if
    end do

end subroutine extract_values

subroutine do_bilinear_interpolation(ncol, nrow, in_array, x, y, zz)
    implicit none

    integer(kind=4), intent(in) :: ncol, nrow
    integer(kind=4) :: i, j
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: in_array(nrow,ncol)
    real(kind=r8), intent(in) :: x, y
    real(kind=r8), intent(out) :: zz
    real(kind=r8) :: fracx, fracy, z11, z21, z12, z22

    j = floor(x)
    i = floor(y)

    fracx = x - j
    fracy = y - i

    call check_i_j(in_array, i, j, nrow, ncol, fracx, fracy, x, y, z11, z21, z12, z22)

    zz = z11 + fracy*(z21-z11) + fracx*(z12-z11) + fracx*fracy*(z11+z22-z21-z12)

end subroutine do_bilinear_interpolation

subroutine check_i_j(in_array, i, j, y_dim, x_dim, fracx, fracy, x, y, z11, z21, z12, z22)
    implicit none

    integer(kind=4), intent(in) :: i, j, y_dim, x_dim
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: in_array(y_dim, x_dim)
    real(kind=r8), intent(in) :: fracx, fracy, x, y
    real(kind=r8) :: z11, z21, z12, z22

    logical :: is_error

    is_error = .false.

    if (j < 1 .or. j > x_dim) then
        print *, 'fuera de rango: j=', j, 'x_dim=', x_dim
        error stop
    end if

    if (i < 1 .or. i > y_dim) then
        print *, 'fuera de rango: i=', i, 'y_dim=', y_dim
        error stop
    end if

    z11 = in_array(i, j)

    ! Round z11

    if (j < x_dim) then
        z12 = in_array(i, j+1)
    end if

    if (i < y_dim) then
        z21 = in_array(i+1, j)
    end if

    if (j < x_dim .and. i < y_dim) then
        z22 = in_array(i+1, j+1)
    end if

    if (x_dim == j) then
        if (fracx == 0) then
            z12 = 0
            z22 = 0
        else
            is_error = .true.
        end if
    end if

    if (y_dim == i) then
        if (fracy == 0) then
            z21 = 0
            z22 = 0
        else
            is_error = .true.
        end if
    end if

    if (j > x_dim .or. i > y_dim) then
        is_error = .true.
    end if

    if (z22 < 0) then
        print *, 'i:', i, 'j:', j
        print *, 'fracx =', fracx
    end if

    if (is_error) then
        print *, 'A source or receiver location is outside the supported landscape'
        print *, '(must be >= 1/2 cell width/length from the edge of the landscape)'
        print *, 'i:', i, 'j:', j, 'x_dim:', x_dim, 'y_dim:', y_dim, 'fracx:', fracx, 'fracy:', fracy, 'x:', x, 'y:', y
        error stop
    end if

end subroutine check_i_j

subroutine calc_vegmax(veg_cut, distance, npts, m2ft, max_veg_loss)
    implicit none

    integer(kind=4), intent(in) :: npts
    integer(kind=4), intent(in) :: veg_cut(npts)
    integer(kind=4) :: count_con, count_heb, count_hwd, i
    integer(kind=4), dimension(1) :: con_val = [2]
    integer(kind=4), dimension(1) :: heb_val = [3]
    integer(kind=4), dimension(1) :: hwd_val = [4]
    integer(kind=4), dimension(1) :: shb_val = [5]
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: distance, m2ft
    real(kind=r8), intent(out) :: max_veg_loss
    real(kind=r8) :: max_con_loss, distance_con, max_hwd_loss, max_heb_loss, distance_hwd

    if (distance < 75) then
        max_veg_loss = 0
    else
        count_con = 0
        do i = 1, npts
            if (any(veg_cut(i) == con_val)) then
                count_con = count_con + 1
            end if
        end do

        if (count_con == 0) then
            max_con_loss = 0
        else
            distance_con = (count_con / real(npts, kind=8)) * (distance / m2ft)
            max_con_loss = 5.2504 * log(distance_con) - 9.8094
            if (max_con_loss < 0.0) then
                max_con_loss = 0.0
            end if
        end if

        count_hwd = 0
        do i = 1, npts
            if (any(veg_cut(i) == hwd_val)) then
                count_hwd = count_hwd + 1
            end if
        end do

        if (count_hwd == 0) then
            max_hwd_loss = 0
        else
            distance_hwd = (count_hwd / real(npts, kind=8)) * (distance / m2ft)
            max_hwd_loss = 6.6224 * log(distance_hwd) - 16.762
            if (max_hwd_loss < 0) then
                max_hwd_loss = 0
            end if
        end if

        count_heb = 0

        do i = 1, npts
            if (any(veg_cut(i) == heb_val) .or. any(veg_cut(i) == shb_val)) then
                count_heb = count_heb + 1
            end if
        end do

        max_heb_loss = 0
        if (count_heb > 0) then
            max_heb_loss = 4
        end if

        max_veg_loss = max_con_loss + max_hwd_loss + max_heb_loss

        if (max_veg_loss > 14) then
            max_veg_loss = 14
        end if
    end if

end subroutine calc_vegmax

subroutine min(values, min_value, n)
    implicit none

    integer(kind=4), intent(in) :: n
    integer(kind=4) :: i
    integer(kind=4), parameter :: r8=selected_real_kind(15,307)

    real(kind=r8), intent(in) :: values(n)
    real(kind=r8) :: min_value

    min_value = values(1)
    do i = 2, n
        if (values(i) < min_value) then
            min_value = values(i)
        end if
    end do

end subroutine min
