function findClean(app)
            % we need to find usable area (rec in shape) in the current scan image
            % in this prototype, it just gives one fixed region into the clean region list. 
            % here is a transplant of the legacy code in exec find clean function. 
            number = app.RegionSize/1e-9;
            
            ele = 0.1;
            picc = (app.ScanIm - min(min(app.ScanIm)));
            [n,m] = size(app.ScanIm);
            bu = edge(picc,'canny',ele);
            cant = bu;
            r = 4;
            for i = 1+r:n-r
                for j = 1+r:m-r
                    if bu(i,j) == 1
                        for ii = -r:r
                            for jj = -r:r
                                i1 = i+ii;
                                j1 = j+jj;
                                if bu(i1,j1) == 1
                                    tt = max([ii,jj]);
                                    kk = min([ii,jj]);
                                    for t = 1:tt
                                        k = round(t*kk/tt);
                                        if tt == ii
                                            cant(i+t,j+k) = 1;
                                        else
                                            cant(i+k,j+t) = 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            fillhole = imfill(cant,'holes');
            dir = fillhole;
            for i = 1:n
                for j = 1:m
                    if cant(i,j) == 1
                        dir(i,j) = 0;
                    end
                end
            end
            
            edgee = fillhole;
            rr = 4;
            for i=1+rr:n-rr
                for j=1+rr:m-rr
                    if dir(i,j)==1
                        edgee(i-rr:i+rr,j-rr:j+rr)=zeros(2*rr+1,2*rr+1);
                    end
                end
            end
            
            r = 2;
            for i = 1+r:n-r
                for j = 1+r:m-r
                    if bu(i,j) == 1
                        for ii = -r:r
                            for jj = -r:r
                                i1 = i+ii;
                                j1 = j+jj;
                                if bu(i1,j1) == 1
                                    tt = max([ii,jj]);
                                    kk = min([ii,jj]);
                                    for t = 1:tt
                                        k = round(t*kk/tt);
                                        if tt == ii
                                            edgee(i+t,j+k) = 1;
                                        else
                                            edgee(i+k,j+t) = 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % clean = zeros(number,number);
            for i = 1:number
                for j = 1:number
                    if (sum(sum(fillhole(round((i-1)*256/number)+1:round(i*256/number),round((j-1)*256/number)+1:round(j*256/number)))) < 0.01*(256/number)^2) && (sum(sum(edgee)) >= 0.001*256^2) && (sum(sum(dir)) <= 0.2*256^2) 
                        % clean(i,j) = 1;
                        position = [(i-0.5)*app.RegionSize/number, (j-0.5)*app.RegionSize/number, 1e-09, 1e-09];
                        cleanRg = regionGenerator(app,position,1,'g');
                        app.CleanList{end+1} = cleanRg;
                    end
                end
            end
        end
