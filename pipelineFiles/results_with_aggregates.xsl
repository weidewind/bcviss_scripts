<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:preserve-space elements="tree" />
<xsl:template match="/*">
  <html>
  <body>
  
  <h3>Your sample contains:<xsl:value-of select="@name"/></h3>
  <xsl:for-each select="tree">
  <table style="font-family:calibri" cellpadding="5">
    <xsl:for-each select="item">
	  <tr>
	    <td><xsl:value-of select="."/></td>
	  </tr>
	</xsl:for-each> 
  </table>	
  </xsl:for-each>
  

  <xsl:for-each select="chromatogram">
  <h3>Chromatogram <xsl:value-of select="@name"/></h3>
   <xsl:text> Length: </xsl:text><xsl:value-of select="@length"/><xsl:text> nt </xsl:text> <br></br>
   <xsl:text> Only sequences with expectation > 0.05 are shown </xsl:text> <br></br>
   <xsl:text> Total expectation of these sequences: </xsl:text> <xsl:value-of select="round(@total_expect * 1000) div 1000"/> <br></br>
   <a href="{@fasta_ref}"><xsl:text> Download</xsl:text> </a><xsl:text>  all sequences in FASTA format</xsl:text> <br></br>
   <br></br>
   		<table style="font-family:calibri" cellpadding="5">
			<xsl:for-each select="tree/item">
				<tr>
					<td><xsl:value-of select="."/></td>
				</tr>
			</xsl:for-each> 
		</table>	
  <xsl:for-each select="candidate">
  <xsl:sort order="descending" select="@expectation"/>
    <h4 style="text-align:left"> 
		Sequence <xsl:value-of select="position()"/> (<xsl:value-of select="@name"/>)
	</h4>
    <a href="{@img_ref}"><xsl:text> View</xsl:text> </a><xsl:text>  phylogenetic tree (SVG format)</xsl:text> <br></br>
	<a href="{@newick_ref}"><xsl:text> Download</xsl:text></a><xsl:text> phylogenetic tree (Newick format)</xsl:text><br></br>
	<br></br>
	<table style="font-family:calibri" border="1" cellpadding="5">
      <tr>
        <th style="text-align:left" bgcolor="#C8D4FA" >Expectation</th>
        <td><xsl:value-of select="round(@expectation * 1000) div 1000"/></td>
      </tr>
      <tr>
        <th style="text-align:left" bgcolor="#C8D4FA" >Length</th>
        <td><xsl:value-of select="@length"/> nt</td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9DB2F5" >Blast best hit taxonomy</th>
        <td><xsl:value-of select="beast/blast/@taxonomy"/></td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9DB2F5" >Blast best hit information</th>
        <td>Identity = <xsl:value-of select="beast/blast/@identity"/> 
		    Coverage = <xsl:value-of select="beast/blast/@coverage"/> 
			E-value = <xsl:value-of select="beast/blast/@e-value"/></td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9DB2F5" >Blast best hit PROKid</th>
        <td><xsl:value-of select="beast/blast/@prokID"/></td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9DB2F5" >Blast best hit PROKname</th>
        <td><xsl:value-of select="beast/blast/@PROKname"/></td>
      </tr>	  
	  <tr>
        <th style="text-align:left" bgcolor="#9acd32" >STAP taxonomy</th>
        <td>
		<table>
		<xsl:for-each select="beast/taxonomy">
		<tr><td><strong><xsl:value-of select="."/></strong></td></tr>
		</xsl:for-each>
		</table>
		</td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9acd32">PROKnames for relatives</th>
        <td><xsl:value-of select="beast/PROKnames"/></td>
      </tr>
	  <tr>
        <th style="text-align:left" bgcolor="#9acd32">STAP confidence</th>
        <td><xsl:value-of select="round(beast/@confidence * 1000) div 1000"/></td>
      </tr>
    </table>
	</xsl:for-each>
	</xsl:for-each>
  </body>
  </html>
</xsl:template>
</xsl:stylesheet>