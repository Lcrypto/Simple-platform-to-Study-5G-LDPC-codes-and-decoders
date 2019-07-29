clear;

BG_sel = 1;

% H_BG = load(['BG',num2str(BG_sel),'_','iLS',num2str(iLS)]);

ldpc_fid = fopen('test.dat','w','n');

fprintf(ldpc_fid, '%8d%8d\n', 2, 8);

%========================================================
for BG_sel=1:2

    H_BG = load(['BG',num2str(BG_sel)]);
    [H_row_num,H_col_num] = size(H_BG);

    H_ABG = zeros(8,H_row_num,H_col_num);

    for i=1:8
        H_ABG(i,:,:) = load(['BG',num2str(BG_sel),'_','iLS',num2str(i)]);
    end

    H_col = zeros(H_col_num,max(sum(H_BG,1))+1);          % Row i contains the indexes of the check nodes in which the variable i is involved
    H_row = zeros(H_row_num,max(sum(H_BG,2))+1);        % Row i contains the indexes of the variable nodes involved in the check i

    for n=1:H_col_num
        tmph_col = find(H_BG(:,n));
        tmp_col_weight = length(tmph_col);
        H_col(n,1) = tmp_col_weight;
        H_col(n,2:tmp_col_weight+1) = tmph_col;
    end

    for m=1:H_row_num
        tmph_row = find(H_BG(m,:));
        tmp_row_weight = length(tmph_row);
        H_row(m,1) = tmp_row_weight;
        H_row(m,2:tmp_row_weight+1) = tmph_row;
    end

    fprintf(ldpc_fid, '%8d%8d\n', H_col_num, H_row_num);
    fprintf(ldpc_fid, '%8d%8d\n', max(sum(H_BG,1)), max(sum(H_BG,2)));
    for m=1:H_row_num
        fprintf(ldpc_fid, '%8d', H_row(m,1));
        for n=1:H_row(m,1)
            fprintf(ldpc_fid, '%8d', H_row(m,n+1)-1);
            for i=1:8
                fprintf(ldpc_fid, '%8d', H_ABG(i,m,H_row(m,n+1)));
            end
        end
        fprintf(ldpc_fid, '\n');
    end
end

% %========================================================
% BG_sel = 2;
% 
% H_BG = load(['BG',num2str(BG_sel)]);
% [H_row_num,H_col_num] = size(H_BG);
% 
% H_ABG1 = zeros(8,H_row_num,H_col_num);
% 
% for i=1:8
%     H_ABG1(i,:,:) = load(['BG',num2str(BG_sel),'_','iLS',num2str(i)]);
% end
% 
% H_col = zeros(H_col_num,max(sum(H_BG,1))+1);          % Row i contains the indexes of the check nodes in which the variable i is involved
% H_row = zeros(H_row_num,max(sum(H_BG,2))+1);        % Row i contains the indexes of the variable nodes involved in the check i
% 
% for n=1:H_col_num
%     tmph_col = find(H_BG(:,n));
%     tmp_col_weight = length(tmph_col);
%     H_col(n,1) = tmp_col_weight;
%     H_col(n,2:tmp_col_weight+1) = tmph_col;
% end
% 
% for m=1:H_row_num
%     tmph_row = find(H_BG(m,:));
%     tmp_row_weight = length(tmph_row);
%     H_row(m,1) = tmp_row_weight;
%     H_row(m,2:tmp_row_weight+1) = tmph_row;
% end
% 
% fprintf(ldpc_fid, '%8d%8d\n', H_col_num, H_row_num);
% for m=1:H_row_num
%     fprintf(ldpc_fid, '%8d', H_row(m,1));
%     for n=1:H_row(m,1)
%         fprintf(ldpc_fid, '%8d', H_row(m,n+1)-1);
%         for i=1:8
%             fprintf(ldpc_fid, '%8d', H_ABG1(i,m,H_row(m,n+1)));
%         end
%     end
% 	fprintf(ldpc_fid, '\n');
% end

fclose(ldpc_fid);