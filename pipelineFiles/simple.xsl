<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/*">
  <html>
	<head>
	<script type="text/javascript" src="C:/Users/weidewind/Documents/Arch/jquery-1.11.1.min.js"></script>
	</head>
  <body>
  <table style="font-family:calibri" cellpadding="5" border="1" cellspacing="0">
  	<tr>
		<th>Chromatogram name</th>
		<th>Expectation</th>
		<th>STAP taxonomy</th>
		<th>PROKMSA names</th>
		<th>Best BLAST hit taxonomy</th>
		<th>Hit identity</th>
		<th>Hit coverage</th>
		<th>Hit PROKMSA name</th>

	</tr>
	<xsl:for-each select="chromatogram">

		<xsl:for-each select="simple">
		<xsl:sort order="descending" select="@expectation"/>
		<tr class="collapsing-toggler">
			<xsl:if test="(count(parent::*/preceding-sibling::*) + 1) mod 2 = 0">
				<xsl:attribute name='BGCOLOR'>#DCEBFA</xsl:attribute> 
			</xsl:if>		
			<xsl:if test="position() = 1">
				<th rowspan="{last()}"><xsl:value-of select="../@name"/></th>
			</xsl:if>
			<td><xsl:value-of select="round(@expectation * 1000) div 1000"/></td>
			<td><xsl:value-of select="@stap_taxonomy"/></td>
			<td><xsl:value-of select="@PROKnames"/></td>
			<td><xsl:value-of select="@blast_taxonomy"/></td>
			<td><xsl:value-of select="@blast_identity"/> </td>
		    <td><xsl:value-of select="@blast_coverage"/> </td>
			<td><xsl:value-of select="@blast_PROKname"/></td>

		</tr>
		
		</xsl:for-each>
		<tr class="collapsing-panel">
			<td colspan = "8">

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
    				    <th style="text-align:left" bgcolor="#9DB2F5" >Blast best hit prokID</th>
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
						<tr><td><xsl:value-of select="."/></td></tr>
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
    				    <td><xsl:value-of select="beast/@confidence"/></td>
    				  </tr>
    				</table>
					</xsl:for-each>
			</td>
		</tr>
	</xsl:for-each>
  </table>
  <script>
  $(document).ready(function(){
	$('.collapsing-panel').hide();
	$('.collapsing-toggler').click(function()
	{
		$(this).nextAll('.collapsing-panel:first').slideToggle(400);
	});
	});
 </script>
  </body>
  </html>
</xsl:template>
</xsl:stylesheet>