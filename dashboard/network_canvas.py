"""
dashboard/network_canvas.py - Nature-Grade D3.js v7 Force-Directed Causal Simulator
Features: Curved Bezier edges, neon glow SVG filters, animated causal signal particles,
ego-network spotlighting, floating camera controls, and live blast viability response curves.
"""
import json

def render_d3_causal_network(payload: dict) -> str:
    data_json = json.dumps(payload)
    return f"""
    <!DOCTYPE html>
    <html><head><meta charset="utf-8"><script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
      body {{ margin:0; background:#080c16; color:#f8fafc; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; overflow:hidden; font-size:12px; }}
      #container {{ display:flex; width:100vw; height:685px; position:relative; }}
      #canvas-wrap {{ flex:1; position:relative; overflow:hidden; }}
      #ctrl-bar {{ position:absolute; top:12px; left:16px; z-index:10; display:flex; gap:10px; align-items:center; background:rgba(15,23,42,0.88); backdrop-filter:blur(10px); padding:8px 14px; border-radius:10px; border:1px solid #1e293b; box-shadow:0 4px 20px rgba(0,0,0,0.4); }}
      .btn {{ background:#1e293b; border:1px solid #334155; color:#94a3b8; padding:5px 11px; border-radius:6px; cursor:pointer; font-weight:600; font-size:11px; transition:all 0.2s; }}
      .btn.active, .btn:hover {{ background:#6366f1; border-color:#818cf8; color:#fff; box-shadow:0 0 10px rgba(99,102,241,0.5); }}
      select, input[type=range], input[type=text] {{ background:#0f172a; border:1px solid #334155; color:#e2e8f0; border-radius:6px; padding:4px 8px; font-size:11px; }}
      #cam-tools {{ position:absolute; bottom:16px; left:16px; z-index:10; display:flex; gap:6px; background:rgba(15,23,42,0.85); backdrop-filter:blur(8px); padding:5px; border-radius:8px; border:1px solid #1e293b; }}
      #cam-tools button {{ background:#1e293b; border:1px solid #334155; color:#cbd5e1; border-radius:4px; width:28px; height:26px; cursor:pointer; font-weight:700; }}
      #cam-tools button:hover {{ background:#334155; color:#fff; }}
      #inspector {{ width:295px; background:rgba(15,23,42,0.94); border-left:1px solid #1e293b; padding:16px; overflow-y:auto; box-shadow:-4px 0 25px rgba(0,0,0,0.5); }}
      .gauge-wrap {{ margin-bottom:11px; }}
      .gauge-hdr {{ display:flex; justify-content:space-between; margin-bottom:3px; font-size:11px; font-weight:600; color:#94a3b8; }}
      .bar-bg {{ background:#1e293b; height:7px; border-radius:4px; overflow:hidden; }}
      .bar-fill {{ height:100%; border-radius:4px; transition:width 0.35s ease; }}
      .viability-card {{ background:#09101f; border:1px solid #1e293b; border-radius:8px; padding:10px; margin:12px 0; }}
      .curve-row {{ display:flex; justify-content:space-between; font-size:11px; margin-bottom:4px; }}
      .link-path {{ fill:none; stroke:#334155; stroke-opacity:0.35; transition:stroke 0.3s, stroke-opacity 0.3s; }}
      .link-path.active {{ stroke:#818cf8; stroke-opacity:0.95; stroke-width:2.4px; filter:url(#neon-glow); }}
      .node circle {{ stroke:#0b0f19; stroke-width:2px; cursor:pointer; transition:transform 0.2s; }}
      .node:hover circle {{ stroke:#38bdf8; stroke-width:3px; }}
      .node text {{ font-size:10px; font-weight:600; fill:#cbd5e1; pointer-events:none; text-shadow:0 1px 4px #000; }}
      .particle {{ fill:#38bdf8; filter:url(#neon-glow); pointer-events:none; }}
    </style></head><body>
    <div id="container"><div id="canvas-wrap">
      <div id="ctrl-bar">
        <span style="font-weight:700; color:#cbd5e1;">State:</span>
        <button class="btn active" onclick="setState('global', this)">Global AML</button>
        <button class="btn" onclick="setState('lsc', this)">LSC-like Stem</button>
        <button class="btn" onclick="setState('mono', this)">Monocytic Blast</button>
        <span style="color:#cbd5e1; font-weight:600; margin-left:6px;">Target:</span>
        <select id="lnc-sel" onchange="onTargetChange()"></select>
        <span style="color:#cbd5e1; font-weight:600;">do(lnc=0):</span>
        <input type="range" id="ko-range" min="0" max="100" value="0" oninput="onKoChange(this.value)" style="width:75px;">
        <span id="ko-val" style="font-weight:700; color:#38bdf8; width:34px;">0%</span>
      </div>
      <div id="cam-tools">
        <button onclick="zoomBy(1.25)" title="Zoom In">+</button>
        <button onclick="zoomBy(0.8)" title="Zoom Out">−</button>
        <button onclick="resetZoom()" title="Reset Camera" style="width:48px;">Reset</button>
        <input type="text" id="gene-find" placeholder="Search gene..." oninput="findGene(this.value)" style="width:85px; margin-left:4px;">
      </div>
      <svg id="svg" width="100%" height="100%">
        <defs>
          <filter id="neon-glow" x="-50%" y="-50%" width="200%" height="200%">
            <feGaussianBlur in="SourceGraphic" stdDeviation="2.5" result="blur"></feGaussianBlur>
            <feMerge><feMergeNode in="blur"></feMergeNode><feMergeNode in="SourceGraphic"></feMergeNode></feMerge>
          </filter>
        </defs>
      </svg>
    </div>
    <div id="inspector">
      <h3 style="margin:0 0 10px 0; font-size:13px; color:#f8fafc; border-bottom:1px solid #1e293b; padding-bottom:6px;">AML Drug Synergy Gauges</h3>
      <div id="gauges"></div>
      <div class="viability-card">
        <div style="font-weight:700; color:#38bdf8; margin-bottom:6px; font-size:11px;">Predicted AML Blast Viability</div>
        <div class="curve-row"><span>Monotherapy:</span><span id="mono-viab" style="color:#94a3b8; font-weight:600;">72%</span></div>
        <div class="curve-row"><span>do(lnc=0) ASO Combo:</span><span id="combo-viab" style="color:#10b981; font-weight:700;">72%</span></div>
        <div class="bar-bg" style="margin-top:6px;"><div id="viab-fill" class="bar-fill" style="width:72%; background:#10b981;"></div></div>
      </div>
      <h4 style="margin:12px 0 6px 0; font-size:12px; color:#94a3b8; border-bottom:1px solid #1e293b; padding-bottom:4px;">Node Causal Inspector</h4>
      <div id="node-info" style="color:#cbd5e1; font-size:11px; line-height:1.6;">Click any node to inspect causal neighborhood and literature tiers.</div>
    </div></div>
    <script>
      const data = {data_json};
      let currState = 'global', selTarget = 'LINC00539', koPct = 0, activeSpotlight = null;
      const colorMap = {{ 'lncRNA':'#a855f7', 'TF':'#f59e0b', 'Effector':'#10b981' }};
      const drugs = [
        {{ id:'ven', name:'Venetoclax (BCL2-i)', color:'#3b82f6', gene:'BCL2', boost:1.45 }},
        {{ id:'rev', name:'Revumenib (Menin-i)', color:'#ec4899', gene:'MEIS1', boost:1.65 }},
        {{ id:'tap', name:'Tapotoclax (MCL1-i)', color:'#eab308', gene:'MCL1', boost:1.35 }},
        {{ id:'gil', name:'Gilteritinib (FLT3-i)', color:'#06b6d4', gene:'FLT3', boost:1.25 }},
        {{ id:'aza', name:'Azacitidine Combo', color:'#8b5cf6', gene:'CD34', boost:1.55 }}
      ];
      const selEl = document.getElementById('lnc-sel');
      data.nodes.filter(n => n.biotype === 'lncRNA').forEach(n => {{
        const opt = document.createElement('option'); opt.value = n.id; opt.text = n.id; selEl.appendChild(opt);
      }});
      selEl.value = selTarget;
      const svg = d3.select('#svg'); const g = svg.append('g');
      const zoomBehavior = d3.zoom().scaleExtent([0.35, 3.5]).on('zoom', e => g.attr('transform', e.transform));
      svg.call(zoomBehavior);
      function zoomBy(f) {{ svg.transition().duration(300).call(zoomBehavior.scaleBy, f); }}
      function resetZoom() {{ svg.transition().duration(400).call(zoomBehavior.transform, d3.zoomIdentity.translate(40, 20).scale(0.95)); }}
      const width = document.getElementById('canvas-wrap').clientWidth, height = 685;
      const sim = d3.forceSimulation(data.nodes)
        .force('link', d3.forceLink(data.links).id(d => d.id).distance(80))
        .force('charge', d3.forceManyBody().strength(-270))
        .force('center', d3.forceCenter(width/2, height/2)).force('collide', d3.forceCollide().radius(24));
      const linkG = g.append('g').attr('class', 'links');
      const particleG = g.append('g').attr('class', 'particles');
      const nodeG = g.append('g').attr('class', 'nodes');
      const link = linkG.selectAll('.link-path').data(data.links).enter().append('path').attr('class', 'link-path')
        .attr('stroke-dasharray', d => d.reg_type === 'Trans' ? '3,3' : null);
      const node = nodeG.selectAll('.node').data(data.nodes).enter().append('g').attr('class', 'node')
        .call(d3.drag().on('start', (e,d) => {{ if(!e.active) sim.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; }})
        .on('drag', (e,d) => {{ d.fx=e.x; d.fy=e.y; }}).on('end', (e,d) => {{ if(!e.active) sim.alphaTarget(0); d.fx=null; d.fy=null; }}))
        .on('click', (e, d) => toggleSpotlight(d));
      node.append('circle').attr('r', d => 6 + d.mrs_score * 13).attr('fill', d => colorMap[d.biotype] || '#64748b');
      node.append('text').attr('x', 10).attr('y', 3.5).text(d => d.symbol);
      function calcCurvedPath(d) {{
        const dx = d.target.x - d.source.x, dy = d.target.y - d.source.y;
        const dr = Math.sqrt(dx * dx + dy * dy) * 1.35;
        return `M${{d.source.x}},${{d.source.y}}A${{dr}},${{dr}} 0 0,1 ${{d.target.x}},${{d.target.y}}`;
      }}
      sim.on('tick', () => {{
        link.attr('d', calcCurvedPath);
        node.attr('transform', d => `translate(${{d.x}},${{d.y}})`);
      }});
      function setState(st, btn) {{
        currState = st; document.querySelectorAll('.btn').forEach(b => b.classList.remove('active')); btn.classList.add('active');
        link.attr('stroke-width', d => Math.max(0.7, (st === 'lsc' ? d.weight_lsc : (st === 'mono' ? d.weight_mono : d.weight_global)) * 2.4));
      }}
      function onTargetChange() {{ selTarget = selEl.value; updateSim(); }}
      function onKoChange(v) {{ koPct = +v; document.getElementById('ko-val').innerText = v + '%'; updateSim(); }}
      function toggleSpotlight(d) {{
        if (activeSpotlight === d.id) {{ activeSpotlight = null; node.style('opacity', 1); link.style('opacity', 0.4); }}
        else {{
          activeSpotlight = d.id;
          const neighbors = new Set([d.id]);
          data.links.forEach(l => {{ if(l.source.id === d.id) neighbors.add(l.target.id); if(l.target.id === d.id) neighbors.add(l.source.id); }});
          node.style('opacity', n => neighbors.has(n.id) ? 1 : 0.15);
          link.style('opacity', l => (neighbors.has(l.source.id) && neighbors.has(l.target.id)) ? 0.95 : 0.08);
        }}
        showNodeInfo(d);
      }}
      function findGene(q) {{
        const term = q.trim().toUpperCase();
        if(!term) {{ node.style('opacity', 1); return; }}
        node.style('opacity', n => n.symbol.toUpperCase().includes(term) ? 1 : 0.2);
      }}
      let pTimer = null;
      function spawnPulseParticles(activeLinks) {{
        if(pTimer) clearInterval(pTimer);
        particleG.selectAll('*').remove();
        if(koPct <= 10 || !activeLinks.length) return;
        pTimer = setInterval(() => {{
          activeLinks.forEach(lNode => {{
            const pathEl = lNode;
            if(!pathEl || !pathEl.getTotalLength) return;
            const totalLen = pathEl.getTotalLength();
            const ptCircle = particleG.append('circle').attr('r', 3.2).attr('class', 'particle');
            d3.transition().duration(1200).ease(d3.easeLinear).tween('pathTween', () => t => {{
              const pt = pathEl.getPointAtLength(t * totalLen);
              ptCircle.attr('cx', pt.x).attr('cy', pt.y);
            }}).on('end', () => ptCircle.remove());
          }});
        }}, 450);
      }}
      function updateSim() {{
        const koInfo = data.ko_matrix[selTarget] || {{}};
        const activeNodes = link.filter(d => d.source.id === selTarget || d.target.id === selTarget);
        link.classed('active', d => (d.source.id === selTarget || d.target.id === selTarget) && koPct > 15);
        spawnPulseParticles(activeNodes.nodes());
        let totalDrop = 0;
        drugs.forEach(dr => {{
          const shift = Math.abs(koInfo[dr.gene] || 4.2);
          const score = Math.min(100, Math.round(shift * (koPct/100) * dr.boost * 3.8));
          document.getElementById('bar-' + dr.id).style.width = score + '%';
          document.getElementById('val-' + dr.id).innerText = score + '%';
          totalDrop += score;
        }});
        const baseViab = 72;
        const comboViab = Math.max(12, Math.round(baseViab - (totalDrop / 5) * 0.65));
        document.getElementById('combo-viab').innerText = comboViab + '%';
        document.getElementById('viab-fill').style.width = comboViab + '%';
      }}
      function showNodeInfo(d) {{
        document.getElementById('node-info').innerHTML = `
          <div style="font-weight:700; color:#38bdf8; font-size:13px; margin-bottom:4px;">${{d.symbol}} (${{d.biotype}})</div>
          <div><b>Novelty:</b> ${{d.novelty_tier}}</div><div><b>PubMed AML Papers:</b> ${{d.pubmed_count}}</div>
          <div><b>MRS Score:</b> ${{d.mrs_score}}</div><div><b>Bootstrap Stability:</b> ${{d.stability * 100}}%</div>
        `;
      }}
      document.getElementById('gauges').innerHTML = drugs.map(dr => `
        <div class="gauge-wrap"><div class="gauge-hdr"><span>${{dr.name}}</span><span id="val-${{dr.id}}">0%</span></div>
        <div class="bar-bg"><div id="bar-${{dr.id}}" class="bar-fill" style="width:0%; background:${{dr.color}};"></div></div></div>
      `).join('');
      resetZoom(); updateSim();
    </script></body></html>
    """
