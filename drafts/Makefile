.PHONY : epsfigs

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

epsfigs : 
	for x in  $$(grep includegr hybrid_zone.tex | sed -e 's/.*{//' | sed -e 's/}//g'); do y=$$(echo $$x | sed -e 's/pdf/eps/'); make $$y; mv $$y for_submission/figs; done
	for x in  $$(grep includegr supplement.tex | sed -e 's/.*{//' | sed -e 's/}//g'); do y=$$(echo $$x | sed -e 's/pdf/eps/'); make $$y; mv $$y for_submission/figs; done
