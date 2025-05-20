function export_pdf_figure(h_num,name,style)

h=figure(h_num);
set(h,'Units','Inches');
pos = get(h,'Position');

if style==1
    set(h,'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]/pos(3)*10);
else
    set(h,'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]/pos(3)*5);
end
set(findobj(gcf,'type','axes'),'FontSize',12)
print(h,name,'-dpdf','-r500','-bestfit')

end