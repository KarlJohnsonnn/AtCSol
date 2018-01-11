%% Funktion zum Auslesen der in Fortran90 produzierten Matrizen
%
%

function [MatrixStruct,m_sparse] = SparseInput(Path)


fileID = fopen(Path,'r');
tline = fgetl(fileID);

% skip headlines
for i=1:6
    tline = fgetl(fileID);
end

% read matrix dimension
tline = fgetl(fileID);
dim = sscanf( tline , '%*19c %d %*3c %d' );
m   = dim(1,1);
n   = dim(2,1);

% read number of nonzero elements
tline = fgetl(fileID);
nonzeros = sscanf( tline , '%*19c %d ' );

% read number of reactions and number of species
tline = fgetl(fileID);
dim = sscanf( tline , '%*19c %d %*3c %d' );
nr  = dim(1);
ns  = dim(2);

% check if CLASSIC or EXTENDED matrix case
tline = fgetl(fileID);

version = sscanf( tline , '%*19c %2c' );
if (strcmp(version,'cl'))
    ncol=ns;
elseif (strcmp(version,'ex'))
    ncol=nr;
end

fclose(fileID);

% allocate arrays
RowInd    = zeros(nonzeros,1);
ColInd    = zeros(nonzeros,1);
Val       = zeros(nonzeros,1);
DiagPtr   = zeros(n,1);
DiagPtr_R = zeros(nr,1);
DiagPtr_C = zeros(ns,1);
ROW_P     = zeros(ns,1);
COL_P     = zeros(ncol,1);
X_P       = 0;
Permu     = zeros(n,2);   % permu and invpermu
LUPermu   = zeros(nonzeros,1);   % permu and invpermu


% another open
fileID = fopen(Path,'r');

holetext = textscan(fileID,'%s');
fclose(fileID);

nCell = length(holetext{1,1}(:));

for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'MATRIX' )) )
        for j=1:nonzeros
            RowInd(j) = str2double(holetext{1,1}{ i+j+cnt   , 1 });
            ColInd(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            Val(j)    = str2double(holetext{1,1}{ i+j+cnt+2 , 1 })+eps;
            cnt = cnt + 2;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'DIAG_PTR' )) )
        for j=1:n
            DiagPtr(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            cnt = cnt + 1;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'DIAG_R_PTR' )) )
        for j=1:nr
            DiagPtr_R(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            cnt = cnt + 1;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'DIAG_C_PTR' )) )
        for j=1:ns
            DiagPtr_C(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            cnt = cnt + 1;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'DIAG_PERMUTATION' )) )
        for j=1:n
            Permu(j,1) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            Permu(j,2) = str2double(holetext{1,1}{ i+j+cnt+2 , 1 });
            cnt = cnt + 2;
        end
        break;
    end
end
% for i=1:nCell
%     cnt = 0;
%     if ( ~isempty(strfind( holetext{1,1}{i,1} , 'LU_PERMUTATION' )) )
%         for j=1:nonzeros
%             LUPermu(j,1) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
%             cnt = cnt + 1;
%         end
%         break;
%     end
% end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'ROW_PTR' )) )
        for j=1:ns
            ROW_P(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            cnt = cnt + 1;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'COL_PTR' )) )
        for j=1:ncol
            COL_P(j) = str2double(holetext{1,1}{ i+j+cnt+1 , 1 });
            cnt = cnt + 1;
        end
        break;
    end
end
for i=1:nCell
    cnt = 0;
    if ( ~isempty(strfind( holetext{1,1}{i,1} , 'X_PTR' )) )
        X_P = str2double(holetext{1,1}{ i+cnt+2 , 1 });
    end
end


% build output stucture
MatrixStruct.path= Path;
MatrixStruct.m   = m;
MatrixStruct.n   = n;
MatrixStruct.nr  = nr;
MatrixStruct.ns  = ns;
MatrixStruct.nnz  = nonzeros;

MatrixStruct.ri  = RowInd;
MatrixStruct.ci  = ColInd;
MatrixStruct.val = Val+eps;

MatrixStruct.dp  = DiagPtr;
MatrixStruct.dpr = DiagPtr_R;
MatrixStruct.dpc = DiagPtr_C;
MatrixStruct.perm = Permu;
MatrixStruct.luperm = LUPermu;
MatrixStruct.rowp = ROW_P;
MatrixStruct.colp = COL_P;
MatrixStruct.xp = X_P;

m_sparse = sparse(RowInd,ColInd,Val+eps,m,n);
end