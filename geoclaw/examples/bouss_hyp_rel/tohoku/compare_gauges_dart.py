from pylab import *
from clawpack.pyclaw.gauges import GaugeSolution


swe_color = 'blue'
HSGN_color = 'green'
SGN_color = 'magenta'
mbHSGN_color =  '#F05800'


event = 'Fujii'
outdir1 = './_output_fujii_SGN'
outdir2 = './_output_fujii_SWE'
outdir3 = './_output_fujii_mbHSGN'
outdir4 = './_output_fujii_HSGN'


fs1 = 18
fs2 = 18

for k,gaugeno in enumerate([21413, 21418]):
    gaugeno1 = gaugeno
    gaugeno2 = gaugeno

    gauge1 = GaugeSolution(gauge_id=gaugeno1, path=outdir1)
    gauge2 = GaugeSolution(gauge_id=gaugeno2, path=outdir2)
    gauge3 = GaugeSolution(gauge_id=gaugeno2, path=outdir3)
    gauge4 = GaugeSolution(gauge_id=gaugeno2, path=outdir4)
    fname = '%s_notide.txt' % gaugeno1
    observed = loadtxt(fname)

    figure(502+k,figsize=(13,4))
    clf()

    plot(gauge2.t/60., gauge2.q[-1,:], color=swe_color, label='SWE',linewidth=2)
    
    plot(gauge3.t/60., gauge3.q[-1,:], color=mbHSGN_color, label='mbHSGN',linewidth=2)
    plot(gauge4.t/60., gauge4.q[-1,:], color=HSGN_color, label='HSGN',linewidth=2, linestyle='--')
    plot(gauge1.t/60., gauge1.q[-1,:], color=SGN_color, label='SGN',linewidth=2,linestyle=':')
    plot(observed[:, 0] / 60., observed[:, 1],
        color='k', linewidth=0.8,
        marker='o', markerfacecolor='none', markeredgecolor='k',
        markeredgewidth=0.8, label='Obs')       
    if gaugeno1 == 21413: 
        xlim(60,180)
        if event == 'Fujii':
            ylim(-0.6,1.2)
        else:
            ylim(-1,1.5)
    if gaugeno1 == 21418: 
        xlim(20,90)
        if event == 'Fujii':
            ylim(-1,2.)
        else:
            ylim(-2,2.5)
    xlabel('Minutes after earthquake',fontsize=fs2)
    legend(loc='upper right',fontsize=fs1, ncol=2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    grid(True)
    ylabel(r'$\eta$ (m)',fontsize=fs2)

    title('Gauge %s comparison using %s source' % (gaugeno1,event),
          fontsize=fs1)
    fname = 'DART%s_%s.pdf' %  (gaugeno1,event)
    savefig(fname, bbox_inches='tight', pad_inches=0)
    print('Created ',fname)
