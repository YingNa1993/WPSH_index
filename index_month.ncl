begin
startyear=1951
endyear=2020
year= ispan(startyear,endyear,1)
nyear=dimsizes(year)
nmonth=12

a=addfile("/data/naying/1_WPSH/hgt.mon.mean.nc","r")
hgtncep=a->hgt 
UNDEF=hgtncep@_FillValue
lat=a->lat
lon=a->lon
nlat=dimsizes(lat)
nlon=dimsizes(lon)

hgt1=(hgtncep(((startyear-1948)*12):((endyear-1948)*12+11),{500},:,:))/10.
hgt=onedtond(ndtooned(hgt1),(/nyear,nmonth,nlat,nlon/))

hgt!0 = "year"
hgt!1 = "month"
hgt!2 = "lat"
hgt!3 = "lon"
hgt&year = year
hgt&lat  = lat
hgt&lon  = lon
hgt&lat@units = "degrees_north"
hgt&lon@units = "degrees_east"

hgt500=hgt(:,:,{10:60},{110:150})
hgt500!2="lat"
hgt500!3="lon"
hgt500&lat  = lat({10:60})
hgt500&lon  = lon({110:150})
hgt500&lat@units = "degrees_north"
hgt500&lon@units = "degrees_east"

nlat1=dimsizes(lat({10:60}))
nlon1=dimsizes(lon({110:150}))

hgt500_1=hgt(:,:,{7.5:62.5},{107.5:152.5})
hgt500_mean9=new((/nyear,nmonth,nlat1,nlon1/),float,UNDEF)
do iy= 0,nyear-1
    do id=  0,nmonth-1 
        do i=1,nlat1
            do j=1,nlon1
                hgt500_mean9(iy,id,i-1,j-1)=avg(hgt500_1(iy,id,(i-1):(i+1),(j-1):(j+1)))
            end do 
        end do 
    end do 
end do 
copy_VarMeta(hgt500, hgt500_mean9)

lon1=(/110,130/)
lon2=(/130,150/)

index_new_GX=new((/2,nyear,nmonth/),float,UNDEF)
index_new_GI=new((/2,nyear,nmonth/),float,UNDEF)


do aa= 0,1
lon11=lon1(aa)
lon22=lon2(aa)
latx=lat({10:60})
lonx=lon({lon11:lon22})
nlatx=dimsizes(latx)
nlonx=dimsizes(lonx)

hgtx=new((/nyear,nmonth,nlatx,nlonx/),float,UNDEF)
hgtx= hgt500(:,:,:,{lon11:lon22})
hgtxm9= hgt500_mean9(:,:,:,{lon11:lon22})
hgtx!2="lat"
hgtx!3="lon"
hgtx&lat  = latx
hgtx&lon  = lonx
hgtx&lat@units = "degrees_north"
hgtx&lon@units = "degrees_east"
copy_VarMeta(hgtx, hgtxm9)


do iy= 0,nyear-1
do id= 0,nmonth-1 
    if (.not.all(ismissing(hgtx(iy,id,:,:)))) then
        
        ;;挑选每一经度上位势高度最大值hgtmax极其所在纬度latmax;;;;;;;;;;;;;;;;;;;;;;;;;;           
        hgtmax   = new(nlonx,float,UNDEF)
        latmax   = new(nlonx,float,UNDEF)
        hgtm9max = new(nlonx,float,UNDEF)
        do ilon=0,nlonx-1
            hgtmax(ilon) = max(hgtx(iy,id,:,ilon))
            nmax = num(hgtx(iy,id,:,ilon).eq.hgtmax(ilon))      

            if (nmax.eq.1) then 
                latmax(ilon) = latx(maxind(hgtx(iy,id,:,ilon)))
                hgtm9max(ilon) = hgtxm9(iy,id,maxind(hgtx(iy,id,:,ilon)),ilon)
            else if (nmax.gt.1) then 
                ;;;;;;;;;测试是否有1个以上的最大值
                ;print(tostring(year(iy))+tostring(id+1)+"  "+tostring(lonx(ilon)))
                ;;;;;;;;;是 使用周围9点平均较大的点
                ilatind := ind(hgtx(iy,id,:,ilon).eq.hgtmax(ilon))
                ;;;;;;;;;测试周围9点平均是否出现相等的情况
                ;hgtmean9max=max(hgtxm9(iy,id,ilatind,ilon))
                ;nmaxmean9=num(hgtxm9(iy,id,ilatind,ilon).eq.hgtmean9max)
                ;if(nmaxmean9.gt.1)then 
                ;print(tostring(year(iy))+tostring(id+1)+"  "+tostring(lonx(ilon)))
                ;end if 
                ;;;;;;;;;;否 未出现过9点平均相等的情况
                ilatindmax=ilatind(maxind(hgtxm9(iy,id,ilatind,ilon)))
                latmax(ilon) = latx(ilatindmax)
                hgtm9max(ilon) = hgtxm9(iy,id,ilatindmax,ilon)
            end if 
            end if 
        end do 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;;挑选{lon11:lon22}范围内较大的nn个hgtmax
        nn= toint(ceil(nlonx/2.))
        ;print("nn="+nn)
        hgtmax_nn =new(nn,float)
        latmax_nn =new(nn,float)
        lonmax_nn =new(nn,float)

        do ii=0,nn-1
            hgtmax_nn(ii) = max(hgtmax)

            nmax=num(hgtmax.eq.hgtmax_nn(ii))
            if (nmax.eq.1) then 
                latmax_nn(ii) = latmax(maxind(hgtmax))
                lonmax_nn(ii) = lonx(maxind(hgtmax))
                hgtmax(maxind(hgtmax)) = UNDEF
            else if (nmax.gt.1) then
                ;;;;;;;;;测试是否有1个以上的最大值
                ;print(tostring(year(iy))+tostring(id+1))
                ;;;;;;;;;是 使用周围9点平均较大的点

                iind := ind(hgtmax.eq.hgtmax_nn(ii))

                ;;;;;;;;;测试周围9点平均是否出现相等的情况
                ;hgtmean9max := max(hgtm9max(iind))
                ;nmaxmean9=num(hgtm9max(iind).eq.hgtmean9max)
                ;if(nmaxmean9.gt.1)then 
                ;print(tostring(year(iy))+tostring(id+1))
                ;end if
                ;;;;;;;;仅出现过1次1972年2月第二第三时相等 不影响结果

                iiind := iind(maxind(hgtm9max(iind)))

                latmax_nn(ii) =latmax( iiind )
                lonmax_nn(ii) =lonx( iiind )
                hgtmax(iiind) = UNDEF

            end if 
            end if 

            

        end do  


        index_new_GX(aa,iy,id)=avg(latmax_nn)
        index_new_GI(aa,iy,id)=avg(hgtmax_nn)
    end if 

end do 
end do 

end do 


GX_west= index_new_GX(0,:,:)
GX_east= index_new_GX(1,:,:)
GI_west= index_new_GI(0,:,:)
GI_east= index_new_GI(1,:,:)

GX=new((/nyear,nmonth/),float)
GI=new((/nyear,nmonth/),float)
GX=(GX_west+GX_east)/2.
GI=(GI_west+GI_east)/2.



;---Open a new NetCDF file to write to
;fout_name = "index_wpsh_new_1951-2020_month.nc"
;system("rm -f " + fout_name)
;fout = addfile(fout_name,"c")
;fout->index_new = GX
;fout->index_new_west = GX_west
;fout->index_new_east = GX_east
;fout->index_new_GI = GI
;fout->index_new_west_GI = GI_west
;fout->index_new_east_GI = GI_east
;delete(fout)


GX=where(ismissing(GX),-9999.,GX)
;GI=where(ismissing(GI),-9999.,GI)
GX_west=where(ismissing(GX_west),-9999.,GX_west)
;GI_west=where(ismissing(GI_west),-9999.,GI_west)
GX_east=where(ismissing(GX_east),-9999.,GX_east)
;GI_east=where(ismissing(GI_east),-9999.,GI_east)

GX1:=new(nyear, string)
;GI1:=new(nyear, string)
     do i = 0, nyear-1
       GX1(i) = str_concat (sprintf ("%9.2f", GX(i,:)))
       ;GI1(i) = str_concat (sprintf ("%9.2f", GI(i,:)))
     end do
asciiwrite ("GX_mo_1951-2020.txt", GX1)
;asciiwrite ("GI_mo_1951-2020.txt", GI1)
GX1:=new(nyear, string)
;GI1:=new(nyear, string)
     do i = 0, nyear-1
       GX1(i) = str_concat (sprintf ("%9.2f", GX_west(i,:)))
       ;GI1(i) = str_concat (sprintf ("%9.2f", GI_west(i,:)))
     end do
asciiwrite ("GX_west_mo_1951-2020.txt", GX1)
;asciiwrite ("GI_west_mo_1951-2020.txt", GI1)
GX1:=new(nyear, string)
;GI1:=new(nyear, string)
     do i = 0, nyear-1
       GX1(i) = str_concat (sprintf ("%9.2f", GX_east(i,:)))
       ;GI1(i) = str_concat (sprintf ("%9.2f", GI_east(i,:)))
     end do
asciiwrite ("GX_east_mo_1951-2020.txt", GX1)
;asciiwrite ("GI_east_mo_1951-2020.txt", GI1)

;;;;txt文件1-70行为1951-2020年，1-12列为月。;;;;;;;



end 